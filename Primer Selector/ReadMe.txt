Primer_Information_2.py采用MPprimer的方法计算引物退火温度。并可以显示最佳匹配
修改：
删除原来已有的rotavirus引物
>Rotavirus_A#F#16
TAATGCTTTTCAGTGGTTGMTGC
>Rotavirus_A#F#23
GGCWTTTAATGCTTTTCAGT
>Rotavirus_A#F#15
MTGCTCAAGATGGAGTCTACWCA

Primers_amplicons_information.py
输入文件为引物序列信息、参考序列信息
输出：
1. 引物信息
2. 引物在参考序列上的位置信息，通过序列比对进行查找，允许错配
3. 引物

PrimerSelector_Amplicon_2.py与PrimerSelector_Amplicon.py类似，但采用Primers_amplicons_information.py中的方法得到引物的位置信息，设定一个阈值(如90%)，引物与序列匹配度大于阈值即可。返回的positon位置为最佳匹配。排除在Amplicon上检索不到的引物.
对PrimerSelector_Amplicon_2.py进行了改写，使用不同的函数将不同功能的模块联系起来，使得结构更加明晰。可以作为今后引物筛选的主要脚本。
可能还需要优化一下速度，目前软件运行速度偏慢。

ProbeSelector.py
输入：
primer_seq, ref_seq, 每个ref_seq对应一组引物序列
处理：
在ref_seq上找到primer_seq的位置，根据一定的规则(长度、Tm值、noG@beginning)等找到候选的探针序列。
输出：
找到的探针序列

用于找到候选引物序列后生成candidate primers

ProbeSelector_2.py
输入：
primer_seq, ref_seq. 每个ref_seq对应多组引物序列
处理：
根据一定条件，找到可能的引物组合，对没一组引物组合，寻找可能存在的候选探针。
输出：
每一对引物的候选探针。

用于比较难设计探针的序列，例如rotavirus的探针，之前找了几组，在特定条件下，都无法生成rotavirus的探针。

PrimerProbeSelector.py
在PrimerSelector_Amplicon_2.py的基础上修改，更加方便。
输入：
primer_seq, ref_seq. 每个ref_seq对应多组引物序列
处理：
根据一定条件，找到可能的引物组合，对没一组引物组合，寻找可能存在的候选探针。
输出：
每一对引物的候选探针。

用于比较难设计探针的序列，例如rotavirus的探针，之前找了几组，在特定条件下，都无法生成rotavirus的探针。

PrimerProbeSelector-primer3.py 在PrimerProbeSelector.py基础上添加了检测hairpin结构的功能，根据primer3 manual的介绍，其计算hairpin的cutoff值为47，暂时不清楚这个47是什么意思。先采用默认的数值。当简并引物中hairpin greater than 25%，则舍弃探针。

 Add hairpinCheck() and homodimerCheck() function to probe, if the  tm of the hairpin structure  is greater than 50, hairpin_number +=1, if hairpin number great than 25% of the total degenerate number, the probe is unqualified.
 If the delta G value of  homodimer structure is less than -9kcal*mol-1, the probe is unqualified

 python PrimerProbeSelector-primer3.py 1primer_tested.fasta 4Refseq_GI_Virus_Plasmids.fasta



python PrimerSelector_Amplicon_2.py GI_Virus_Primers_0727.fasta Refseq_GI_Virus_Plasmids.fasta
python PrimerSelector_Amplicon_2.py chosen_GI_Virus_Primers_0727.fasta Refseq_GI_Virus_Plasmids.fasta
python PrimerSelector_Amplicon_2.py mergedPrimers2.fasta Refseq_GI_Virus_Plasmids.fasta
python PrimerSelector_Amplicon_2.py mergedPrimers3.fasta Refseq_GI_Virus_Plasmids.fasta
python PrimerSelector_Amplicon_2.py GI_Virus_4_primers.fasta Refseq_GI_Virus_Plasmids.fasta



python PrimerSelector_Amplicon_2.py test_primers.fasta test_input_region.fasta
python ProbeSelector.py Primers_Ordered_Ref.fasta Refseq_GI_Virus_Plasmids.fasta
python ProbeSelector.py test_probe.fasta Refseq_GI_Virus_Plasmids.fasta

#目前只能用于单对引物，用于多对引物会出错
python PrimerProbeSelector.py Rotavirus_Primers.fasta Rotavirus_Refseq.fasta

python PrimerProbeSelector.py Norovirus_g2.fasta Norovirus_g2_ref.fasta.fasta

修改PrimerProbeSelector.py，可以用于多对引物。
PrimerProbeSelector-primer3_2.py在PrimerProbeSelector-primer3.py的基础上修改，且可用于多对引物。先测试一下。
python PrimerProbeSelector-primer3.py GI_5_multiplex.fasta 5Refseq_GI_Virus_Plasmids.fasta
python PrimerProbeSelector-primer3_2.py GI_5_multiplex.fasta 5Refseq_GI_Virus_Plasmids.fasta
python PrimerProbeSelector.py GI_5_multiplex.fasta 5Refseq_GI_Virus_Plasmids.fasta

改写PrimerProbeSelector.py，修改两个方面：
1. 可以生成的探针区域包括正向序列和反向互补序列
2. 探针的命名以ref-seq@f/r@strat-index@end-index来命名
python PrimerProbeSelector.py updates_primers_e9_candi.fasta 6Refseq_GI_Virus_Plasmids_updates.fasta

修改之后将其整合到PrimerProbeSelector-primer3_2中.py
python PrimerProbeSelector_primer3_2.py GI_4_multiplex.fasta 4Refseq_GI_Virus_Plasmids_updates.fasta

python PrimerProbeSelector_primer3_2.py updates_primers_e9_candi.fasta 6Refseq_GI_Virus_Plasmids_updates.fasta

PrimerSelector_Amplicon_RT.py 用于使用realtime PCR对引物进行筛选
PrimerSelector_Amplicon_CE.py 用于使用毛细管电泳法对引物进行筛选