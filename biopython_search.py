import scispacy
import spacy
import pandas as pd
import streamlit as st
from Bio import Entrez, Medline
from io import BytesIO
from scispacy.linking import EntityLinker

# Load the SciSpacy model and UMLS linker (optional)
nlp = spacy.load("en_core_sci_sm")  # You can use en_core_sci_lg for a larger, more accurate model
linker = EntityLinker(resolve_abbreviations=True, name="umls")

# Function to extract entities from abstracts using SciSpacy
def extract_entities_scispacy(text):
    doc = nlp(text)
    chemicals = []
    mesh_terms = []
    
    for ent in doc.ents:
        # Filter relevant entities like chemicals, diseases, genes, etc.
        if ent.label_ in ["CHEMICAL", "DISEASE", "ENTITY", "GENE_OR_GENOME"]:  # You can refine this based on your needs
            chemicals.append(ent.text)
        
        # Linking entities to UMLS or other ontologies for MeSH terms
        if linker is not None:
            for umls_ent in ent._.umls_ents:
                linked_entity = linker.kb.cui_to_entity[umls_ent[0]]
                if "MeSH" in linked_entity.definition:  # Example of identifying MeSH terms
                    mesh_terms.append(linked_entity.canonical_name)
    
    return mesh_terms, chemicals

# Function to fetch articles and extract entities
def fetch_and_process_articles(query, num_articles, email):
    Entrez.email = email
    st.write(f"Searching PubMed with query: {query}")
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=num_articles)
        result = Entrez.read(handle)
        handle.close()

        ids = result['IdList']
        if not ids:
            st.write("No articles found matching the criteria.")
            return []

        handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        articles = list(records)
        handle.close()

        # Extract entities from abstracts
        processed_data = []
        for article in articles:
            title = article.get('TI', 'No title available')
            abstract = article.get('AB', 'No abstract available')
            
            # Process text with SciSpacy
            mesh_terms, chemicals = extract_entities_scispacy(abstract)
            
            processed_data.append({
                'Title': title,
                'Abstract': abstract,
                'MeSH Terms': mesh_terms,
                'Chemicals': chemicals
            })
        
        return pd.DataFrame(processed_data)
    
    except Exception as e:
        st.write(f"An error occurred: {e}")
        return []

# Function to generate Excel file in memory
def save_to_excel(df):
    output = BytesIO()
    
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False)
    
    output.seek(0)  # Move the buffer position to the beginning
    return output

# Streamlit UI for user inputs
st.title("PubMed Research Navigator with SciSpacy")
email = st.text_input("Enter your email (for PubMed access):")
search_term = st.text_input("Enter the general search term:")
num_articles = st.number_input("Enter the number of articles to fetch:", min_value=1, max_value=1000, value=10)

if st.button("Fetch and Analyze Articles"):
    if email and search_term:
        query = f"({search_term})"
        articles_df = fetch_and_process_articles(query, num_articles, email)
        
        if not articles_df.empty:
            st.write("Extracted Data:")
            st.write(articles_df)
            
            # Allow download of processed data
            excel_data = save_to_excel(articles_df)
            st.download_button(
                label="Download Extracted Data",
                data=excel_data.getvalue(),
                file_name="extracted_pubmed_data.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        else:
            st.write("No articles fetched.")
    else:
        st.write("Please enter valid input.")

# Display copyright information at the bottom of the app
st.write("""
    ### Copyright Information
    Copyright (c) 2024 Anurag Pande
    
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
    to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
    IN THE SOFTWARE.
""")
