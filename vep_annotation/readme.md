# 将VEP注释结果转换为BGICG注释结果

## 整体框架拆分为以下几个步骤：

1）从原始VEP中提取需要信息，从vcf文件中注释出来变异支持数目；

2）各类库文件注释

3）其他各种基于上述的条件判断

4）其他


#
### 从VEP中提取的列包括
1）Gene ： EGFR

2）cHGVS : c.451C>A

3）pHGVS ：p.P151T

4）pHGVS2 ：p.PRO151=

5）ExIn_ID ：EX5；

6）vep_Function ：vep原始注释结果

7）Function：vep简化后的结果

8）vep2bgicg_Function: vep简化后的结果对应到bgicg结果上

9）SIFT:

10）polyphen2：

11）chr

12）start: 使用标准格式

13）end:

14）ref：

15）alt：

16）Muttye: 根据Ref和Alt列的碱基数据输出

17）Genotype: G/T 看下格式，是否需要反向互补之类的

18）Transcript：NM_000546.5

19）Protein：NP_000537.3

20）Strand：正负链 +-

21）Flank侧翼序列：看下怎么来的，写脚本获取

22）rsID：直接提取

23）BI_MutType： 蛋白注释功能、（来自vep的VARIANT_CLASS）

24）clinvar：直接提取

25）
Local_AF
Local_Hom_AF
1000G_AF
1000G_EAS_AF
G1000_AMR_AF
G1000_AFR_AF
G1000_EUR_AF
G1000_SAS_AF
ESP6500
ExAC_AF
ExAC_EAS_AF
ExAC_EAS_Hom_AF
EXAC_AFR_AF
ExAC_AFR_Hom_AF
EXAC_AMR_AF
ExAC_AMR_Hom_AF
EXAC_FIN_AF
ExAC_FIN_Hom_AF
EXAC_NFE_AF
ExAC_NFE_Hom_AF
EXAC_SAS_AF
ExAC_SAS_Hom_AF
