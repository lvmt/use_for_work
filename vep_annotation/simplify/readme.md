# 此模块的功能为将VEP注释结果转化为BGI解读需要格式 


## 关于Func信息转换的逻辑

### 一、vep功能简化
1） vep原始Func列为单一字段，则vep_simple_func直接继承

2） vep原始Func列为2个子Func，进行以下判断：

    2.1） 如果这2个子Func属于不同的基因区域（exon，intro，utr等），则vep_simple_func为span类型

    2.2） 如果这2个子Func属于相同的基因区域，则根据Func优先级信息，选取优先级较高的func

3） vep原始Func列大于2个子Func，则属于span类型
    
#
**<font color=green>bug**

bug修复：2021年12月31日：

vep原始Func列大于2个子Func，不一定属于span类型

示例：chr15:66679686delA, 

func字段信息：frameshift_variant,start_lost,start_retained_variant</font>
    
#

**<font color=red>基于上述bug的发现，将调整整个vep功能简化逻辑</font>**
  - 一： vep 原始Func列为单一字段，则vep_simple_func直接继承
  - 二： 若vep原始列大于1个，则进行以下判断：
     
    2.1） 如果这2个子Func属于不同的基因区域（exon，intro，utr等），则vep_simple_func为span类型

    2.2） 如果这2个子Func属于相同的基因区域，则根据Func优先级信息，选取优先级较高的func
    
#

### 二、vep简化功能转化为BGI解读功能
1） 争对vep种特殊的功能进行额外处理

    1.1） protein_altering_variant
    1.2） coding_sequence_variant

2）判断为span的结果中，存在部分特例

    2.1）争对span的类型，在进行vep2bgi之前，会先进行span的矫正

3）除去1、2两种情况外，直接根据vep2bgi对应关系进行转换



## Function转换中存在的重大问题记录
1） 争对protein_altering_variant需要进行特殊处理

2） 争对coding_sequence_variant需要进行特殊处理 

争对protein_altering_variant的处理，详细见函数handle_protein_altering_variant_from_hgvsp;
处理逻辑由lvmengting和白帅帅共同商讨，由于目前涉及到的案例较少，部分不确定情况返回null

coding_sequence_variant单独出现的案例，目前也比较少，暂时争对出现的案例做了逻辑判断，其余未出现情况返回空值



#

