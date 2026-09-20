import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, MeltingTemp as mt
import pandas as pd

# ضبط إعدادات الصفحة
st.set_page_config(page_title="BioMedical Sequence & Diagnostic Analyzer", page_icon="🧬", layout="wide")

st.title("🧬 BioMedical Sequence & Diagnostic Analyzer (v2.0)")
st.write("An engineering-focused bioinformatics tool bridging clinical sequence analysis with diagnostic parameters.")

# تقسيم التطبيق إلى 3 أجزاء رئيسية
tab1, tab2, tab3 = st.tabs([
    "📊 Basic Analysis & Translation", 
    "🔬 Diagnostic Primer Design (PCR)", 
    "🩺 Clinical Mutation Detector"
])

# ---------------------------------------------------------
# TAB 1: التحليل الأساسي والترجمة
# ---------------------------------------------------------
with tab1:
    st.header("📊 Sequence Analysis & Central Dogma")
    sequence_input = st.text_area("Enter DNA sequence:", "ATGCGATCGATCGATCGATCGATCGATCGATCTAG", height=120, key="t1_seq")
    cleaned_seq = "".join(sequence_input.split()).upper()

    if cleaned_seq:
        if set(cleaned_seq).issubset(set("ATCGN")):
            dna_seq = Seq(cleaned_seq)
            col1, col2 = st.columns(2)
            col1.metric("Sequence Length (bp)", len(dna_seq))
            col2.metric("GC-Content", f"{gc_fraction(dna_seq)*100:.2f}%")
            
            st.subheader("Transcription & Translation")
            st.code(f"RNA: {dna_seq.transcribe()}", language="text")
            protein = dna_seq.translate()
            st.code(f"Protein: {protein}", language="text")
            
            # رسم توزيع الأحماض الأمينية
            aa_counts = {aa: protein.count(aa) for aa in set(protein) if aa != "*"}
            if aa_counts:
                df = pd.DataFrame(list(aa_counts.items()), columns=['Amino Acid', 'Frequency']).set_index('Amino Acid')
                st.bar_chart(df)
        else:
            st.error("⚠️ Invalid DNA sequence! Please use A, T, C, G bases only.")

# ---------------------------------------------------------
# TAB 2: تصميم وتطوير بادئات الـ PCR (هندسة أجهزة التشخيص)
# ---------------------------------------------------------
with tab2:
    st.header("🔬 Diagnostic PCR Primer Thermodynamics")
    st.write("Analyze primer suitability for thermal cyclers and diagnostic qPCR devices.")
    
    primer_input = st.text_input("Enter Diagnostic Primer Sequence (18-30 bp):", "ATGCGATCGATCGATCGATC", key="t2_primer")
    cleaned_primer = "".join(primer_input.split()).upper()
    
    if cleaned_primer:
        if set(cleaned_primer).issubset(set("ATCG")):
            primer_seq = Seq(cleaned_primer)
            
            # حساب درجة حرارة الانصهار بالتنبوء الحراري (Thermodynamical Nearest-Neighbor)
            tm_val = mt.Tm_NN(primer_seq)
            gc_val = gc_fraction(primer_seq) * 100
            
            c1, c2, c3 = st.columns(3)
            c1.metric("Primer Length", f"{len(primer_seq)} bp")
            c2.metric("GC Ratio", f"{gc_val:.1f}%")
            c3.metric("Melting Temp ($T_m$)", f"{tm_val:.2f} °C")
            
            # تقييم ملاءمة البادئ لأجهزة الفحص الحراري
            st.subheader("Diagnostic Suitability Verdict:")
            if 18 <= len(primer_seq) <= 30 and 55 <= tm_val <= 65 and 40 <= gc_val <= 60:
                st.success("✅ **Optimal Primer:** Excellent parameters for standard diagnostic PCR assays.")
            else:
                st.warning("⚠️ **Sub-optimal Primer:** Parameters fall outside ideal qPCR device limits ($T_m$: 55-65°C, GC: 40-60%, Length: 18-30 bp).")
        else:
            st.error("⚠️ Invalid primer bases! Use A, T, C, G only.")

# ---------------------------------------------------------
# TAB 3: كاشف الطفرات السريرية (Clinical Mutation Detector)
# ---------------------------------------------------------
with tab3:
    st.header("🩺 Clinical Mutation & Variant Detector")
    st.write("Compare a patient sequence against a reference sequence to detect single nucleotide mutations (SNPs).")
    
    ref_seq = st.text_input("Reference DNA Sequence (Normal):", "ATGCGATCGATCGATC", key="ref_seq").strip().upper()
    pat_seq = st.text_input("Patient DNA Sequence (Sample):", "ATGCGACCGATCGATC", key="pat_seq").strip().upper()
    
    if ref_seq and pat_seq:
        if len(ref_seq) != len(pat_seq):
            st.warning("⚠️ Sequences must be of equal length for direct alignment comparison.")
        else:
            mutations = []
            for i in range(len(ref_seq)):
                if ref_seq[i] != pat_seq[i]:
                    mutations.append((i + 1, ref_seq[i], pat_seq[i]))
            
            if mutations:
                st.error(f"🚨 Detected {len(mutations)} mutation(s) in the patient sample:")
                mut_df = pd.DataFrame(mutations, columns=["Position (bp)", "Reference Base", "Patient Base"])
                st.table(mut_df)
            else:
                st.success("✅ No point mutations detected! Patient sequence matches reference 100%.")
