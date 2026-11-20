$g="C:\Users\Administrator\WorkBuddy\2026-11-01-11-02-06\glaucoma_R_pipeline"
$dst="$g\data\raw\lincs\Cmap_QCeddata"
$names = @(
"adipose_ASC.C_ph2_QCed.txt.gz","adipose_ASC_ph1_QCed.txt.gz","adipose_ASC_ph2_QCed.txt.gz",
"breast_BT20_ph1_QCed.txt.gz","breast_BT20_ph2_QCed.txt.gz","breast_HME1_ph2_QCed.txt.gz",
"breast_HS578T_ph1_QCed.txt.gz","breast_HS578T_ph2_QCed.txt.gz","breast_MCF10A_ph1_QCed.txt.gz",
"breast_MCF10A_ph2_QCed.txt.gz","breast_MCF7_ph1_QCed.txt.gz","breast_MCF7_ph2_QCed.txt.gz",
"breast_MDAMB231_ph1_QCed.txt.gz","breast_MDAMB231_ph2_QCed.txt.gz","breast_SKBR3_ph1_QCed.txt.gz",
"breast_SKBR3_ph2_QCed.txt.gz","cns_MNEU.E_ph2_QCed.txt.gz","cns_NEU_ph1_QCed.txt.gz",
"cns_NEU_ph2_QCed.txt.gz","cns_NPC.CAS9_ph2_QCed.txt.gz","cns_NPC.TAK_ph2_QCed.txt.gz",
"cns_NPC_ph1_QCed.txt.gz","cns_NPC_ph2_QCed.txt.gz","digestive_AGS_ph1_QCed.txt.gz",
"digestive_CL34_ph1_QCed.txt.gz","digestive_HCT116_ph1_QCed.txt.gz","digestive_HT115_ph1_QCed.txt.gz",
"digestive_HT29_ph1_QCed.txt.gz","digestive_HT29_ph2_QCed.txt.gz","digestive_LOVO_ph1_QCed.txt.gz",
"digestive_MDST8_ph1_QCed.txt.gz","digestive_NCIH508_ph1_QCed.txt.gz","digestive_RKO_ph1_QCed.txt.gz",
"digestive_SNU1040_ph1_QCed.txt.gz","digestive_SNUC4_ph1_QCed.txt.gz","digestive_SNUC5_ph1_QCed.txt.gz",
"digestive_SW480_ph1_QCed.txt.gz","digestive_SW620_ph1_QCed.txt.gz","digestive_SW948_ph1_QCed.txt.gz",
"hematopoietic_CD34_ph1_QCed.txt.gz","hematopoietic_CD34_ph2_QCed.txt.gz","hematopoietic_HL60_ph1_QCed.txt.gz",
"hematopoietic_HS27A_n14_QCed.txt.gz","hematopoietic_JURKAT_ph1_QCed.txt.gz","hematopoietic_JURKAT_ph2_QCed.txt.gz",
"hematopoietic_NOMO1_ph1_QCed.txt.gz","hematopoietic_PL21_ph1_QCed.txt.gz","hematopoietic_SKM1_ph1_QCed.txt.gz",
"hematopoietic_THP1_ph1_QCed.txt.gz","hematopoietic_U266_ph1_QCed.txt.gz","hematopoietic_U937_ph1_QCed.txt.gz",
"hematopoietic_WSUDLCL2_ph1_QCed.txt.gz","liver_HEPG2_ph1_QCed.txt.gz","liver_HEPG2_ph2_QCed.txt.gz",
"liver_HUH7_ph1_QCed.txt.gz","liver_PHH_ph1_QCed.txt.gz","lung_A549_ph1_QCed.txt.gz",
"lung_A549_ph2_QCed.txt.gz","lung_CORL23_ph1_QCed.txt.gz","lung_DV90_ph1_QCed.txt.gz",
"lung_H1299_ph1_QCed.txt.gz","lung_HCC15_ph1_QCed.txt.gz","lung_HCC515_ph1_QCed.txt.gz",
"lung_HCC515_ph2_QCed.txt.gz","lung_NCIH1694_ph1_QCed.txt.gz","lung_NCIH1836_ph1_QCed.txt.gz",
"lung_NCIH2073_n361_QCed.txt.gz","lung_NCIH596_ph1_QCed.txt.gz","lung_SKLU1_ph1_QCed.txt.gz",
"lung_T3M10_ph1_QCed.txt.gz","musculoskeletal_A673_ph1_QCed.txt.gz","musculoskeletal_SKB_ph1_QCed.txt.gz",
"musculoskeletal_SKL.C_ph2_QCed.txt.gz","musculoskeletal_SKL_ph2_QCed.txt.gz","ovary_COV644_ph1_QCed.txt.gz",
"ovary_EFO27_ph1_QCed.txt.gz","ovary_OV7_ph1_QCed.txt.gz","ovary_RMGI_ph1_QCed.txt.gz",
"ovary_RMUGS_ph1_QCed.txt.gz","ovary_TYKNU_ph1_QCed.txt.gz","pancreas_YAPC_ph2_QCed.txt.gz",
"prostate_LNCAP_ph2_QCed.txt.gz","prostate_PC3_ph1_QCed.txt.gz","prostate_PC3_ph2_QCed.txt.gz",
"prostate_VCAP_ph1_QCed.txt.gz","skin_A375_ph1_QCed.txt.gz","skin_A375_ph2_QCed.txt.gz",
"skin_MCH58_ph1_QCed.txt.gz","skin_SKMEL1_ph1_QCed.txt.gz","skin_SKMEL28_ph1_QCed.txt.gz",
"uterus_HEC108_ph1_QCed.txt.gz","uterus_HELA_ph2_QCed.txt.gz","uterus_JHUEM2_ph1_QCed.txt.gz",
"uterus_SNGM_ph1_QCed.txt.gz"
)
$todo = @($names | Where-Object { -not (Test-Path (Join-Path $dst $_)) -or ((Get-Item (Join-Path $dst $_) -ErrorAction SilentlyContinue).Length -lt 5000) })
Write-Output "todo: $($todo.Count)"
$batchSize = 8
for ($i = 0; $i -lt $todo.Count; $i += $batchSize) {
  $batch = $todo[$i..([Math]::Min($i+$batchSize-1, $todo.Count-1))]
  $jobs = foreach ($n in $batch) {
    $u = "https://media.githubusercontent.com/media/konumat/Trans-Phar/master/Cmap_QCeddata/$n"
    $out = Join-Path $dst $n
    Start-Job -ScriptBlock { param($u,$out) & curl.exe -sL --fail --retry 6 --retry-delay 3 -o $out $u } -ArgumentList $u,$out
  }
  $jobs | Wait-Job -Timeout 900 | Out-Null
  $jobs | Remove-Job -Force -ErrorAction SilentlyContinue
  $done = (Get-ChildItem $dst -Filter "*.gz" | Where-Object { $_.Length -gt 5000 }).Count
  Write-Output ("progress: {0}/94" -f $done)
}
$files = Get-ChildItem $dst -Filter "*.gz"
Write-Output ("final: {0} files, {1:N1} MB" -f $files.Count, (($files | Measure-Object Length -Sum).Sum/1MB))
