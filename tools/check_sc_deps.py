#!/usr/bin/env python
# monocle3 / CellChat 硬依赖 vs 已装包 对照
import re
import subprocess

BASE = r"C:/Users/Administrator/WorkBuddy/2026-11-01-11-02-06/glaucoma_R_pipeline/tools"
RSCRIPT = r"G:/R/R-4.5.2/bin/Rscript.exe"

out = subprocess.run([RSCRIPT, "--slave", "-e", "cat(rownames(installed.packages()), sep='\\n')"],
                     capture_output=True, text=True)
installed = set(out.stdout.split())

for pkg_file, name in [("monocle3_src/DESCRIPTION", "monocle3"), ("CellChat_src/DESCRIPTION", "CellChat")]:
    desc = open(f"{BASE}/{pkg_file}", encoding="utf-8").read()
    desc = re.sub(r"\n\s+", " ", desc)
    print(f"===== {name} =====")
    for field in ["Depends", "Imports"]:
        m = re.search(field + r":\s*(.*?)\n[A-Z]\w*:", desc + "\nZ:")
        if not m:
            continue
        deps = [re.sub(r"\s*\(.*\)", "", d).strip() for d in m.group(1).split(",")]
        deps = [d for d in deps if d and d != "R"]
        missing = [d for d in deps if d not in installed]
        print(f"{field}: {len(deps)} 项，缺 {len(missing)}: {', '.join(missing) if missing else '无'}")
    print()
