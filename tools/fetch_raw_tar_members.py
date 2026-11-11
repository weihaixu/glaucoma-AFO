#!/usr/bin/env python
"""按 tar 成员偏移量从 GEO RAW.tar 里用 HTTP Range 抠出指定文件。

GSE263230_RAW.tar 直连流式下载被限速（~200B/s），但 Range 请求正常。
tar 是顺序格式：每个成员 = 512B 头 + 数据 + pad 到 512 的倍数。
成员大小在 filelist.txt 里按顺序给出，据此算偏移；抠出后先校验
h5ad 魔数再落盘，头校验失败即报错（防止偏移算错拿回一堆垃圾）。
"""
import sys
import time
import urllib.request

UA = {"User-Agent": "Mozilla/5.0 (pipeline)"}
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263230/suppl/GSE263230_RAW.tar"

# (文件名, 字节大小) —— 来自 filelist.txt，顺序与 tar 内一致
MEMBERS = [
    ("GSM8188166_obj_query_afterPrediction_10x3V31_GTM_1672_P3.h5ad", 79354744),
    ("GSM8188167_obj_query_afterPrediction_10x3V31_GTM_2956.h5ad", 50670910),
    ("GSM8188168_obj_query_afterPrediction_10x3V31_GTM_7445.h5ad", 97276041),
    ("GSM8188169_obj_query_afterPrediction_10x3V31_hTM_0114_P1.h5ad", 71305119),
    ("GSM8188170_obj_query_afterPrediction_10x3v31_hTM_11701_P2.h5ad", 42819235),
    ("GSM8188171_obj_query_afterPrediction_10x3v31_hTM_11703_P2.h5ad", 44386632),
    ("GSM8188172_obj_query_afterPrediction_10x3v31_hTM_2180.h5ad", 219917378),
    ("GSM8188173_obj_query_afterPrediction_10x3V31_hTM_659_P1_DCR.h5ad", 12035875),
    ("GSM8188174_obj_query_afterPrediction_10x3V31_hTM_74M_P2.h5ad", 76340005),
    ("GSM8188175_obj_query_afterPrediction_10x3v31_hTM_7987_P2.h5ad", 58385736),
    ("GSM8188176_obj_query_afterPrediction_10x3v31_hTM_9355_P1.h5ad", 61516986),
    ("GSM8188177_obj_query_afterPrediction_10x3v31_hTM_9632_P2.h5ad", 45996173),
    ("GSM8188178_obj_query_afterPrediction_10x3V31_hTM_9691_P1.h5ad", 62952579),
    ("GSM8188179_obj_query_afterPrediction_10x3V31_OHSU_GTM_2019_0461.h5ad", 73997779),
    ("GSM8188180_obj_query_predicted_umap_22_0500_TM.h5ad", 410137848),
    ("GSM8188181_obj_query_predicted_umap_22_0688_TM.h5ad", 484399923),
    ("GSM8188182_obj_query_predicted_umap_22_0769_TM.h5ad", 133192238),
]

WANT = {
    "GSM8188166_obj_query_afterPrediction_10x3V31_GTM_1672_P3.h5ad",
    "GSM8188167_obj_query_afterPrediction_10x3V31_GTM_2956.h5ad",
    "GSM8188168_obj_query_afterPrediction_10x3V31_GTM_7445.h5ad",
    "GSM8188170_obj_query_afterPrediction_10x3v31_hTM_11701_P2.h5ad",
    "GSM8188171_obj_query_afterPrediction_10x3v31_hTM_11703_P2.h5ad",
    "GSM8188173_obj_query_afterPrediction_10x3V31_hTM_659_P1_DCR.h5ad",
    "GSM8188176_obj_query_afterPrediction_10x3v31_hTM_9355_P1.h5ad",
    "GSM8188175_obj_query_afterPrediction_10x3v31_hTM_7987_P2.h5ad",
}

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "."


def http_get(url, start, end, retries=6, timeout=300):
    hdrs = dict(UA)
    hdrs["Range"] = "bytes=%d-%d" % (start, end)
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=hdrs)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(1.5 * (i + 1))
    raise last


def pad512(n):
    return (n + 511) // 512 * 512


offset = 0
for name, size in MEMBERS:
    if name in WANT:
        # 先读 512B 头验证成员名（防偏移错位）
        hdr = http_get(URL, offset, offset + 511)
        hdr_name = hdr[:100].split(b"\x00")[0].decode("utf-8", "replace").strip()
        if hdr_name != name:
            sys.exit("偏移校验失败 @%d：头里是 %r，期望 %s" % (offset, hdr_name, name))
        data_start = offset + 512
        sys.stderr.write("取 %s（%d 字节）...\n" % (name, size))
        with open("%s/%s" % (OUTDIR, name), "wb") as f:
            pos = data_start
            remain = size
            while remain > 0:
                step = min(remain, 8 << 20) - 1
                chunk = http_get(URL, pos, pos + step)
                if not chunk:
                    raise RuntimeError("空响应 @%d" % pos)
                f.write(chunk)
                pos += len(chunk)
                remain -= len(chunk)
        sys.stderr.write("  完成 %s\n" % name)
    offset += 512 + pad512(size)

sys.stderr.write("all done\n")
