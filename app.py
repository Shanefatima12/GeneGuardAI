import streamlit as st
import pandas as pd
from Bio.Data.CodonTable import standard_dna_table
from datetime import datetime
from difflib import SequenceMatcher
import io

# ============================================================
# GENEGUARD AI | Educational genetic variant explorer
# ============================================================
st.set_page_config(page_title="GeneGuard AI", page_icon="🧬", layout="wide")

st.markdown("""
<style>
.stApp {background:linear-gradient(135deg,#071421 0%,#0b2540 55%,#102d4b 100%);color:#eaf4ff}
[data-testid="stSidebar"] {background:#07111d}
h1,h2,h3 {color:#70e5ff!important}
[data-testid="stMetric"] {background:#153653;padding:12px;border-radius:12px;border:1px solid #28618a}
div[data-testid="stAlert"] {border-radius:12px}
.base {display:inline-block;min-width:24px;text-align:center;padding:5px 3px;margin:2px;border-radius:5px;font-family:monospace;font-weight:bold;color:#061321}
.A {background:#66e3b4}.T {background:#ffb86b}.G {background:#7bd5ff}.C {background:#c6a5ff}
.changed {outline:3px solid #ffdf57;box-shadow:0 0 8px #ffdf57}
.smallnote {color:#b9d4e8;font-size:0.9rem}
</style>
""", unsafe_allow_html=True)

CODONS = standard_dna_table.forward_table
STOPS = set(standard_dna_table.stop_codons)
VALID = set("ATGC")
AA_NAMES = {
    "A":"Alanine","R":"Arginine","N":"Asparagine","D":"Aspartic acid",
    "C":"Cysteine","Q":"Glutamine","E":"Glutamic acid","G":"Glycine",
    "H":"Histidine","I":"Isoleucine","L":"Leucine","K":"Lysine",
    "M":"Methionine","F":"Phenylalanine","P":"Proline","S":"Serine",
    "T":"Threonine","W":"Tryptophan","Y":"Tyrosine","V":"Valine","*":"Stop"
}

def clean(seq):
    return "".join(str(seq).upper().split()).replace("U","T")

def is_valid(seq):
    return bool(seq) and set(seq) <= VALID

def translate(seq):
    """Translate coding DNA from its first base; incomplete trailing codon is ignored."""
    seq = seq[:len(seq)//3*3]
    return "".join(CODONS.get(seq[i:i+3], "*") for i in range(0,len(seq),3))

def codon_list(seq):
    return [seq[i:i+3] for i in range(0,len(seq)-2,3)]

def mutation_events(ref, var):
    """SequenceMatcher-based event list; positions refer to 1-based reference coordinates."""
    matcher = SequenceMatcher(a=ref, b=var, autojunk=False)
    events = []
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            continue
        rpart, vpart = ref[i1:i2], var[j1:j2]
        if tag == "replace":
            if len(rpart) == len(vpart):
                for offset, (rb, vb) in enumerate(zip(rpart, vpart)):
                    events.append({"Type":"Substitution","Reference position (1-based)":i1+offset+1,
                                  "Reference":rb,"Variant":vb})
            else:
                events.append({"Type":"Complex replacement","Reference position (1-based)":i1+1,
                              "Reference":rpart or "—","Variant":vpart or "—"})
        elif tag == "delete":
            events.append({"Type":"Deletion","Reference position (1-based)":i1+1,
                          "Reference":rpart,"Variant":"—"})
        elif tag == "insert":
            events.append({"Type":"Insertion","Reference position (1-based)":i1+1,
                          "Reference":"—","Variant":vpart})
    return events

def classify_variant(ref, var, ref_protein, var_protein):
    delta = len(var)-len(ref)
    if delta:
        if delta % 3:
            return "Frameshift (candidate)", "Length change is not divisible by 3; this can shift downstream codon grouping. Exact biological interpretation requires alignment and transcript context."
        return "In-frame insertion/deletion", "Length change is divisible by 3, so the reading frame may be preserved; the inserted/deleted amino acids can still affect protein function."
    if ref_protein == var_protein:
        return "Synonymous / no protein change", "The sequences differ but the translated amino-acid strings are identical in this reading frame."
    if "*" in var_protein and "*" not in ref_protein:
        return "Nonsense / premature stop", "The variant translation contains a stop where the reference did not. Confirm the position and transcript before interpreting."
    if "*" in ref_protein and "*" in var_protein and ref_protein.index("*") != var_protein.index("*"):
        return "Stop-codon difference", "The position of a stop codon differs between translations."
    return "Missense / protein-altering", "At least one translated amino acid differs. This alone does NOT establish that a variant is harmful or pathogenic."

def colored_sequence(seq, other, label):
    html = []
    for i, base in enumerate(seq):
        differs = i >= len(other) or base != other[i]
        html.append(f'<span class="base {base} {"changed" if differs else ""}" title="{label} position {i+1}">{base}</span>')
    return "".join(html)

# Sidebar navigation
st.sidebar.title("🧬 GeneGuard AI")
page = st.sidebar.radio("Navigate", ["Variant Analyzer","Disease Explorer","Learning Hub","About & Sources"])
st.sidebar.info("Educational prototype • Not for clinical diagnosis")
st.sidebar.caption("DNA Day Project | Courage • Care • Code")

st.title("🧬 GeneGuard AI")
st.subheader("Interactive Genetic Variant & Disease Awareness Platform")
st.write("Explore how coding-DNA changes can alter a protein, inspect mutation events, and learn about documented disease-associated examples.")
st.divider()

if page == "Variant Analyzer":
    st.header("01 / DNA Variant Analyzer")
    st.caption("Enter aligned coding-DNA sequences (5′→3′) that begin in the same reading frame. This simplified tool does not perform genomic alignment.")
    example = st.selectbox("Load an educational example", [
        "Sickle-cell HBB (missense)", "Synonymous example", "Nonsense example",
        "Frameshift example", "Custom sequences"
    ])
    examples = {
        "Sickle-cell HBB (missense)":("GAG","GTG"),
        "Synonymous example":("GAA","GAG"),
        "Nonsense example":("TGG","TGA"),
        "Frameshift example":("ATGGAGCTG","ATGGAAGCTG"),
        "Custom sequences":("ATGGAGCTG","ATGGAGCTG")
    }
    default_ref, default_var = examples[example]
    col1, col2 = st.columns(2)
    with col1:
        ref_input = st.text_area("Reference coding DNA (5′→3′)", default_ref, height=105, key="ref_input")
    with col2:
        var_input = st.text_area("Variant coding DNA (5′→3′)", default_var, height=105, key="var_input")
    analyze = st.button("🔬 Analyze sequence variant", type="primary", use_container_width=True)

    if analyze:
        ref, var = clean(ref_input), clean(var_input)
        if not is_valid(ref) or not is_valid(var):
            st.error("Please enter both sequences using only A, T, G and C (spaces/newlines are okay).")
            st.session_state.pop("analysis_result", None)
        elif len(ref) < 3 or len(var) < 3:
            st.error("Each sequence must contain at least one complete codon (3 bases).")
            st.session_state.pop("analysis_result", None)
        else:
            st.session_state.analysis_result = (ref, var)

    if "analysis_result" in st.session_state:
        ref, var = st.session_state.analysis_result
        events = mutation_events(ref, var)
        ref_protein, var_protein = translate(ref), translate(var)
        category, explanation = classify_variant(ref, var, ref_protein, var_protein)

        m1,m2,m3,m4 = st.columns(4)
        m1.metric("Reference bases",len(ref)); m2.metric("Variant bases",len(var))
        m3.metric("Detected event groups",len(events)); m4.metric("Length change",f"{len(var)-len(ref):+}")
        st.subheader("Mutation map")
        st.markdown("**Reference**")
        st.markdown(colored_sequence(ref,var,"Reference"),unsafe_allow_html=True)
        st.markdown("**Variant**")
        st.markdown(colored_sequence(var,ref,"Variant"),unsafe_allow_html=True)
        st.caption("A/T/G/C are color-coded • Gold outline marks positional differences (insertions/deletions can shift positions)")
        st.subheader("Variant event table")
        st.dataframe(pd.DataFrame(events) if events else pd.DataFrame(columns=["Type","Reference position (1-based)","Reference","Variant"]),hide_index=True,use_container_width=True)

        st.subheader("DNA → mRNA → Protein")
        c1,c2 = st.columns(2)
        with c1:
            st.markdown("**Reference**")
            st.code(f"Coding DNA: {ref}\n mRNA:       {ref.replace('T','U')}\n Protein:    {ref_protein or '(empty)'}")
        with c2:
            st.markdown("**Variant**")
            st.code(f"Coding DNA: {var}\n mRNA:       {var.replace('T','U')}\n Protein:    {var_protein or '(empty)'}")

        st.subheader("Variant impact classification")
        st.success(category)
        st.write(explanation)
        st.warning("Sequence-based educational classification only. It is not a pathogenicity score or clinical interpretation.")

        st.subheader("Interactive codon explorer")
        ref_codons, var_codons = codon_list(ref), codon_list(var)
        common_count = min(len(ref_codons),len(var_codons))
        if common_count:
            index = st.selectbox("Choose codon position",list(range(common_count)),
                format_func=lambda i:f"Codon {i+1}: {ref_codons[i]} → {var_codons[i]}")
            rc, vc = ref_codons[index], var_codons[index]
            ra, va = CODONS.get(rc,"Stop"), CODONS.get(vc,"Stop")
            st.info(f"Codon {index+1} | Reference: {rc} → {AA_NAMES.get(ra,ra)} ({ra}) | Variant: {vc} → {AA_NAMES.get(va,va)} ({va})")
        else:
            st.info("No shared complete codon positions to display.")

        st.subheader("Amino-acid differences")
        aa_rows=[]
        for i in range(max(len(ref_protein),len(var_protein))):
            ra = ref_protein[i] if i<len(ref_protein) else "—"
            va = var_protein[i] if i<len(var_protein) else "—"
            if ra != va:
                aa_rows.append({"AA position (approx.)":i+1,"Reference":f"{AA_NAMES.get(ra,ra)} ({ra})","Variant":f"{AA_NAMES.get(va,va)} ({va})"})
        st.dataframe(pd.DataFrame(aa_rows) if aa_rows else pd.DataFrame(columns=["AA position (approx.)","Reference","Variant"]),hide_index=True,use_container_width=True)

        report = f"""GENEGUARD AI — EDUCATIONAL VARIANT REPORT
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}
Reference coding DNA (5'→3'): {ref}
Variant coding DNA (5'→3'): {var}
Reference mRNA: {ref.replace('T','U')}
Variant mRNA: {var.replace('T','U')}
Reference protein: {ref_protein}
Variant protein: {var_protein}
Classification: {category}
Explanation: {explanation}
Event groups: {events}
Amino-acid differences: {aa_rows}

DISCLAIMER: Educational sequence comparison only. Not a clinical diagnosis,
validated variant interpretation, pathogenicity prediction, or genetic counseling.
"""
        st.download_button("⬇ Download educational report (.txt)",report,"GeneGuardAI_variant_report.txt","text/plain")
        st.caption("Tip for demo: use the sickle-cell example and explain the coding change GAG → GTG and Glu (E) → Val (V).")

elif page == "Disease Explorer":
    st.header("02 / Disease Explorer")
    disease_data = {
        "Sickle cell disease":{"gene":"HBB","protein":"β-globin","variant":"Classic coding example: GAG → GTG (traditional mature β-globin numbering: Glu6Val; HGVS precursor numbering commonly p.Glu7Val)","impact":"Missense substitution: glutamic acid (E) → valine (V)","summary":"Hemoglobin S can polymerize under deoxygenated conditions, contributing to red-cell sickling. Clinical severity depends on genotype and other factors.","ref":"https://www.ncbi.nlm.nih.gov/books/NBK1377/"},
        "Beta-thalassemia":{"gene":"HBB","protein":"β-globin","variant":"Many substitutions, splice variants, and insertions/deletions","impact":"Variant-specific; can reduce or abolish β-globin production","summary":"Reduced β-globin synthesis can cause anemia and related complications. This is a group of disorders with diverse molecular causes.","ref":"https://medlineplus.gov/genetics/condition/beta-thalassemia/"},
        "Cystic fibrosis":{"gene":"CFTR","protein":"CFTR chloride/bicarbonate channel","variant":"Many CFTR variants; examples include deletion of phenylalanine at position 508 (variant-specific nomenclature)","impact":"Variant-specific protein processing/function effects","summary":"CFTR dysfunction alters ion transport and affects multiple organs, especially lungs and digestive system.","ref":"https://medlineplus.gov/genetics/gene/cftr/"},
        "Phenylketonuria (PKU)":{"gene":"PAH","protein":"Phenylalanine hydroxylase","variant":"Many PAH variants","impact":"Variant-specific enzyme activity changes","summary":"Reduced phenylalanine hydroxylase activity can increase phenylalanine; early detection and management are important.","ref":"https://medlineplus.gov/genetics/gene/pah/"}
    }
    disease = st.selectbox("Choose a documented educational example",list(disease_data))
    d = disease_data[disease]
    st.markdown(f"### {disease}")
    a,b = st.columns(2)
    a.metric("Gene",d["gene"]); b.metric("Protein",d["protein"])
    st.markdown("**Variant/example**"); st.write(d["variant"])
    st.markdown("**Molecular consequence**"); st.write(d["impact"])
    st.markdown("**Disease mechanism (overview)**"); st.write(d["summary"])
    st.markdown(f"**Reference:** [{d['ref']}]({d['ref']})")
    st.info("A disease association does not mean every variant in that gene causes disease. Interpretation depends on the exact variant, zygosity, transcript, evidence, and clinical context.")

elif page == "Learning Hub":
    st.header("03 / Genetics Learning Hub")
    st.markdown("### Central dogma")
    st.success("DNA → RNA → Protein")
    st.write("Transcription produces RNA from DNA information. Translation reads mRNA codons to build a polypeptide.")
    st.markdown("### Mutation types")
    st.markdown("""
- **Substitution:** one nucleotide is replaced by another.
- **Synonymous:** codon changes but encoded amino acid stays the same (in the selected frame).
- **Missense:** substitution changes an amino acid.
- **Nonsense:** a substitution introduces a stop codon.
- **Insertion/deletion:** bases are added or removed; a non-multiple-of-three length change can shift the reading frame.
- **In-frame indel:** length change divisible by three; may still affect protein function.
""")
    st.markdown("### Quick knowledge check")
    question = st.radio("Which type of substitution changes an encoded amino acid?",["Synonymous","Missense","Neither","All substitutions"])
    if st.button("Check answer"):
        if question == "Missense": st.success("Correct! A missense variant changes the encoded amino acid.")
        else: st.info("Try again: missense is the term for a substitution that changes an amino acid.")
    st.markdown("### Mini activity")
    codon = st.selectbox("Explore a coding DNA codon",["GAG","GTG","GAA","TGG","TGA","TTT"])
    aa = CODONS.get(codon,"Stop")
    st.info(f"{codon} → {codon.replace('T','U')} (mRNA) → {AA_NAMES.get(aa,aa)} ({aa})")

else:
    st.header("04 / About & Sources")
    st.write("GeneGuard AI is a student-built educational prototype for coding-DNA comparison, translation, simplified variant classification, and genetics awareness.")
    st.markdown("""
### What it can do
- Compare two user-provided coding-DNA strings
- Display sequence differences and event groups
- Translate coding DNA in a user-selected assumed reading frame (currently starts at base 1)
- Compare translated amino-acid strings
- Provide educational disease examples and a downloadable text report

### Important limitations
- No genomic alignment or reference genome coordinates
- No transcript selection, splice prediction, or validated clinical variant calling
- Simplified indel handling; complex variants can be represented imperfectly
- No pathogenicity score, diagnosis, treatment recommendation, or genetic counseling
- Protein consequences depend on correct sequence, strand, transcript, and reading frame

### Scientific references
- [MedlinePlus Genetics — HBB](https://medlineplus.gov/genetics/gene/hbb/)
- [GeneReviews — Sickle Cell Disease](https://www.ncbi.nlm.nih.gov/books/NBK1377/)
- [NCBI — Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- [MedlinePlus Genetics — CFTR](https://medlineplus.gov/genetics/gene/cftr/)
- [MedlinePlus Genetics — PAH](https://medlineplus.gov/genetics/gene/pah/)
- [MedlinePlus Genetics — Beta-thalassemia](https://medlineplus.gov/genetics/condition/beta-thalassemia/)
""")
st.sidebar.caption("GeneGuard AI • Student educational project")
