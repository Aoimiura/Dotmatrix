"""
23261084 三浦葵衣
"""

import streamlit as st
#10
from Bio import SeqIO

#ファイルからの読み込み
def dotmatrix(f1,f2,win):
   record1 = next(SeqIO.parse(f1,"fasta"))
   record2 = next(SeqIO.parse(f2,"fasta"))

#配列取り出し
   seq1 = record1.seq
   seq2 = record2.seq

#スライス　：　10～19までの10塩基を取り出し
   print(seq1[10:20])
#結果　：　AAATATCCAG

#逆相補鎖
   print(seq1[10:20].reverse_complement())
#結果　：　CTGGATATTT

#翻訳　アミノ酸列への変換
   print(seq1[10:19].translate())
#結果　：　YTF

   import numpy as np
   import matplotlib.pyplot as plt
   #win = 10
   image=np.zeros((500,500))
   image=np.zeros(((len(seq2)-win+1),(len(seq1)-win+1))) #500*500の零行列


   len1=len(seq1)-win+1
   len2=len(seq2)-win+1

   width=500
   height=500

   image=np.zeros((height,width))

#ハッシュ
   hash = {}
#ハッシュに含まれていない？
   subseq1 = "AGATG" #部分配列
   x=102 #部分配列の位置
   for x in range(len1):
     sub1 = seq1[x:x+win]
     if sub1 not in hash:
      hash[sub1]=[]
     hash[sub1].append(x) #位置xを追加


#ハッシュに含まれている？
   subseq2="AGATG"
   y=102
   for y in range(len2):
    sub2 = seq2[y:y+win]
    py=int(y/len2*height)
   if sub2 in hash:
    for x in hash[sub2]:
        px=int(x/len1*width)
        image[py,px]=1
    #subseq2の位置とhash[subseq2]のそれぞれの位置がドット

   plt.imshow(image,extent=(1,len1,len2,1),cmap="Grays")
#plt.show

   st.pyplot(plt)
   
st.title("Dot matrix")
file1=st.sidebar.file_uploader("Sequence file 1:")
file2=st.sidebar.file_uploader("Sequence file 2:")
win = st.sidebar.slider("Window size:",4,100,10)

from io import StringIO

if file1 and file2:
    with StringIO(file1.getvalue().debode("utf-8")) as f1,\
        StringIO(file2.getvalue().decode("utf-8")) as f2:
       dotmatrix(f1,f2,win)
