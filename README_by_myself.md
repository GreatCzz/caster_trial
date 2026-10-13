# CASTER_TRIAL

## Overwiew

CASTER_TRIAL（以下简称 TRIAL），是基于四叶树得分的位点法系统发育树构建工具，基于CASTER工具进行开发。  
此工具旨在减少同源性错误对CASTER的影响，主要有三点改动：  
1、仅计算包含参考物种的四叶树得分。  
2、引入多参考基因组。  
3、基于与参考基因组距离对四叶树得分进行加权（wTRIAL）。  

### 四叶树得分算法

### 引入多参考基因组

### 基于参考基因组距离对四叶树得分加权

## 安装/编译

## 使用样例

'''  
/example/align_ref_a.fas  
>a  
AGCT  
>b  
GGCT  
>c  
GACT  
>d  
AGCG  
>e  
AGCC  
  
/example/align_ref_b.fas  
>a  
CCAA  
>b  
CCTA  
>c  
CCTA  
>d  
CGAA  
>e  
CGAA  
  
/input.txt  
/example/align_ref_a.fas a  
/example/align_ref_b.fas b  
  
TRIAL -i input.txt -o result.tre  
'''