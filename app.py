import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import pandas as pd

st.set_page_config(page_title="DNA Sequence Analyzer", page_icon="🧬")

st.title("🧬 DNA & Protein Sequence Analyzer")
st.write("أداة سريعة لتحليل السلاسل الجينية، حساب نسبة GC، وترجمة الحمض النووي إلى بروتين.")

# ادخال النص
sequence_input = st.text_area("أدخل سلسلة DNA هنا (مثال: ATGC...):", "ATGCGATCGATCGATCGATCGATCGATCG", height=150)

# تنظيف السلسلة
cleaned_seq = "".join(sequence_input.split()).upper()

if cleaned_seq:
    # التحقق من صحة القواعد
    valid_bases = set("ATCGN")
    if set(cleaned_seq).issubset(valid_bases):
        dna_seq = Seq(cleaned_seq)
        
        st.write("---")
        st.subheader("📊 1. التحليل الأساسي")
        col1, col2 = st.columns(2)
        col1.metric("طول السلسلة (bp)", len(dna_seq))
        col2.metric("نسبة GC-Content", f"{gc_fraction(dna_seq)*100:.2f}%")
        
        st.subheader("🔄 2. النسخ والترجمة (Transcription & Translation)")
        rna_seq = dna_seq.transcribe()
        protein_seq = dna_seq.translate()
        
        st.write("**سلسلة RNA:**")
        st.code(str(rna_seq), language="text")
        
        st.write("**سلسلة البروتين (Amino Acids):**")
        st.code(str(protein_seq), language="text")
        
        st.subheader("📈 3. توزيع الأحماض الأمينية")
        aa_counts = {aa: protein_seq.count(aa) for aa in set(protein_seq) if aa != "*"}
        if aa_counts:
            df = pd.DataFrame(list(aa_counts.items()), columns=['الحامض الأميني', 'التكرار']).set_index('الحامض الأميني')
            st.bar_chart(df)
    else:
        st.error("⚠️ السلسلة تحتوي على حروف غير صحيحة! يرجى إدخال (A, T, C, G) فقط.")
