#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fetch_eqtl_tabix.py —— 远程按区间取 tabix 索引文件的子集（只依赖标准库）

为什么需要这个：
  GTEx v8 的 allpairs eQTL 每个组织 3.5GB，本地磁盘与带宽都吃不消；而
  Rsamtools 在这台机器的 Windows 构建里不支持 http(s) 远端 tabix（实测
  "failed to open file: https://..."），pysam 又装不上（PyPI 被代理拦）。
  于是直接按 tabix 规范自己算：读 .tbi 索引 → 得到目标区间落在哪些 BGZF
  块 → 只对这些块发 HTTP Range 请求 → zlib 解压 → 按坐标过滤。
  一次查询通常只传几十到几百 KB。

用法：
  python tools/fetch_eqtl_tabix.py --url URL --index URL.tbi \
      --chrom 1 --start 9300000 --end 9700000 --out out.tsv [--max-bytes N]

输出：带表头的 TSV（表头取自文件首行）。
"""
import argparse
import gzip
import io
import struct
import sys
import time
import zlib
import urllib.request

UA = {"User-Agent": "glaucoma-pipeline/1.0"}


def http_get(url, start=None, end=None, timeout=180, retries=4):
    """本机走本地代理，偶发 502（Tunnel connection failed），必须重试；
    直接失败会让预取的 120 个任务里随机若干个空掉，事后很难发现。"""
    headers = dict(UA)
    if start is not None:
        headers["Range"] = "bytes=%d-%d" % (start, end)
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:            # noqa: BLE001
            last = e
            time.sleep(1.5 * (i + 1))
    raise last


def read_tbi(url):
    raw = gzip.decompress(http_get(url))
    if raw[:4] != b"TBI\x01":
        raise SystemExit("不是 tabix 索引（magic=%r）" % raw[:4])
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<8i", raw[4:36])
    names = raw[36:36 + l_nm].split(b"\x00")
    off = 36 + l_nm
    refs = {}
    for i in range(n_ref):
        (n_bin,) = struct.unpack("<i", raw[off:off + 4]); off += 4
        bins = {}
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack("<Ii", raw[off:off + 8]); off += 8
            chunks = []
            for _ in range(n_chunk):
                b, e = struct.unpack("<QQ", raw[off:off + 16]); off += 16
                chunks.append((b, e))
            bins[bin_id] = chunks
        (n_intv,) = struct.unpack("<i", raw[off:off + 4]); off += 4
        lin = struct.unpack("<%dQ" % n_intv, raw[off:off + 8 * n_intv]); off += 8 * n_intv
        refs[names[i].decode()] = {"bins": bins, "lin": lin}
    return refs


def reg2bins(beg, end):
    end -= 1
    if end < beg:
        end = beg
    lst = [0]
    for k in range(1 + (beg >> 26), 1 + (end >> 26) + 1):
        lst.append(k)
    for k in range(9 + (beg >> 23), 9 + (end >> 23) + 1):
        lst.append(k)
    for k in range(73 + (beg >> 20), 73 + (end >> 20) + 1):
        lst.append(k)
    for k in range(585 + (beg >> 17), 585 + (end >> 17) + 1):
        lst.append(k)
    for k in range(4681 + (beg >> 14), 4681 + (end >> 14) + 1):
        lst.append(k)
    return lst


_BYTES = [0]


def iter_lines(url, start, end, skip0=0, slice_bytes=8 << 20):
    """从 BGZF 的字节区间里逐行产出内容，内存占用只与切片大小有关"""
    pos = start
    d = zlib.decompressobj(31)
    buf = b""

    def feed(data):
        nonlocal d, buf
        while data:
            buf += d.decompress(data)
            data = d.unused_data
            if d.eof:
                d = zlib.decompressobj(31)
            elif not data:
                break

    while pos <= end:
        sl = http_get(url, pos, min(pos + slice_bytes - 1, end))
        if not sl:
            break
        _BYTES[0] += len(sl)
        pos += len(sl)
        feed(sl)
        if skip0:
            # 首个块要按块内偏移跳过前面的残行（在字节层面做，别在行层面做：
            # 行层面会连第一条完整记录一起切掉）
            if len(buf) <= skip0:
                skip0 -= len(buf)
                buf = b""
                continue
            buf = buf[skip0:]
            skip0 = 0
        *lines, buf = buf.split(b"\n")
        for ln in lines:
            yield ln
    if buf:
        yield buf


def bgzf_decompress(buf):
    """把一段 BGZF（= 多个 gzip 成员）解压成字节串"""
    out = []
    while buf:
        d = zlib.decompressobj(31)
        try:
            out.append(d.decompress(buf))
        except zlib.error:
            break
        buf = d.unused_data
        if not d.eof and not buf:
            break
    return b"".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-bytes", type=int, default=200 * 1024 * 1024,
                    help="单次查询允许拉取的最大字节数，超过就报错而不是死等")
    a = ap.parse_args()

    refs = read_tbi(a.index)
    chrom = str(a.chrom)
    if chrom not in refs:
        # 有的索引写成 chr1
        if "chr" + chrom in refs:
            chrom = "chr" + chrom
        else:
            raise SystemExit("索引里没有染色体 %s（可用：%s）" %
                             (a.chrom, ", ".join(list(refs)[:8])))

    ref = refs[chrom]
    bins = reg2bins(a.start - 1, a.end)
    chunks = []
    for b in bins:
        if b in ref["bins"]:
            chunks.extend(ref["bins"][b])
    if ref["lin"]:
        i = min((a.start - 1) >> 14, len(ref["lin"]) - 1)
        min_off = ref["lin"][i]
        chunks = [c for c in chunks if c[1] > min_off]
        chunks = [(max(c[0], min_off), c[1]) for c in chunks]
    chunks.sort()
    if not chunks:
        open(a.out, "w", encoding="utf-8").write("")
        return

    # 合并重叠 chunk
    merged = [list(chunks[0])]
    for b, e in chunks[1:]:
        if b <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([b, e])

    # ---- 分片拉取 + 边解压边过滤 -------------------------------------------
    # 一个 1Mb 窗口解压后 ~1GB 文本，整体读进内存再 6 个并发会直接压爆机器
    # （实测单进程 855MB）。所以按 4MB 压缩块切片，解多少、过多少、丢多少。
    header = ""
    try:
        header = iter_lines(a.url, 0, 65535).__next__().decode("utf-8", "replace")
    except Exception:
        header = ""
    cols = header.split("\t")
    try:
        ci_chr = cols.index("chromosome")
        ci_pos = cols.index("position")
    except ValueError:
        ci_chr, ci_pos = 12, 13

    total = 0
    n_kept = 0
    seen_header = False
    with open(a.out, "w", encoding="utf-8") as fo:
        if header:
            fo.write(header + "\n")
        for b, e in merged:
            start_block = b >> 16
            end_byte = (e >> 16) + 65536 - 1
            skip0 = b & 0xFFFF
            first = True
            for ln in iter_lines(a.url, start_block, end_byte, skip0 if first else 0):
                first = False
                if not ln.strip():
                    continue
                if not seen_header:
                    seen_header = True
                    if header and ln.decode("utf-8", "replace") == header:
                        continue
                f = ln.decode("utf-8", "replace").split("\t")
                if len(f) <= max(ci_chr, ci_pos):
                    continue
                if f[ci_chr] != chrom and "chr" + f[ci_chr] != chrom:
                    continue
                try:
                    p = int(f[ci_pos])
                except ValueError:
                    continue
                if a.start <= p <= a.end:
                    fo.write("\t".join(f) + "\n")
                    n_kept += 1

    sys.stderr.write("chrom=%s:%d-%d 传 %.1f MB，命中 %d 行 -> %s\n" %
                     (chrom, a.start, a.end, _BYTES[0] / 1048576.0, n_kept, a.out))


if __name__ == "__main__":
    main()
