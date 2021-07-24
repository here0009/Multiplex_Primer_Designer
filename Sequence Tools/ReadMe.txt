SequenceMerge.py用于合并序列，
输入：需要合并的fasta文件
处理：将所有独特的fasta序列储存在一个dictionary中，如果序列与已有序列相同，则排除掉。如果序列名称与已有序列相同，则改变名称，新名称为原有名称后面加一个number
输出：会打印出对哪些序列做了修改。输出合并后的文件。
python SequenceMerge.py .\
SequenceSplit.py用于分离序列
将一个含有多条fasta序列的文件分离为多个fasta文件，每个含有一条fasta序列。