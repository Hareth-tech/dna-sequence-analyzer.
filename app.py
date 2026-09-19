import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import pandas as pd

st.set_page_config(page_title="DNA Sequence Analyzer", page_icon="🧬")

st.title("🧬 DNA & Protein Sequence Analyzer")
st.write("A fast tool for DNA sequence analysis, GC-content calculation, and protein translation.")

# Input text area
sequence_input = st.text_area("Enter DNA sequence here (e.g., ATGC...):", "ATGCGATCGATCGATCGATCGATCGATCGATCTAG", height=150)

# Clean sequence
cleaned_seq = "".join(sequence_input.split()).upper()

if cleaned_seq:
    # Validate nucleotides
    valid_bases = set("ATCGN")
    if set(cleaned_seq).issubset(valid_bases):
        dna_seq = Seq(cleaned_seq)
        
        st.write("---")
        st.subheader("📊 1. Basic Analysis")
        col1, col2 = st.columns(2)
        col1.metric("Sequence Length (bp)", len(dna_seq))
        col2.metric("GC-Content", f"{gc_fraction(dna_seq)*100:.2f}%")
        
        st.subheader("🔄 2. Transcription & Translation")
        rna_seq = dna_seq.transcribe()
        protein_seq = dna_seq.translate()
        
        st.write("**RNA Sequence:**")
        st.code(str(rna_seq), language="text")
        
        st.write("**Protein Sequence (Amino Acids):**")
        st.code(str(protein_seq), language="text")
        
        st.subheader("📈 3. Amino Acid Distribution")
        aa_counts = {aa: protein_seq.count(aa) for aa in set(protein_seq) if aa != "*"}
        if aa_counts:
            df = pd.DataFrame(list(aa_counts.items()), columns=['Amino Acid', 'Frequency']).set_index('Amino Acid')
            st.bar_chart(df)
    else:
        st.error("⚠️ Invalid sequence! Please enter valid DNA bases (A, T, C, G) only.")
