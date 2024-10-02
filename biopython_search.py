import os
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from io import BytesIO
from transformers import AutoTokenizer, AutoModelForTokenClassification, pipeline
from Bio import Entrez, Medline
from collections import Counter

# Load BioBERT Tokenizer and Model for NER
tokenizer = AutoTokenizer.from_pretrained("dmis-lab/biobert-base-cased-v1.1")
model = AutoModelForTokenClassification.from_pretrained("dmis-lab/biobert-base-cased-v1.1")
ner_pipeline = pipeline("ner", model=model, tokenizer=tokenizer)

# Function to perform NER on biomedical text using BioBERT
def bio_ner(text):
    # Ensure the text is not empty
    if not text.strip():  # Check for empty or whitespace-only text
        return []

    # Truncate the text to 512 tokens to avoid tensor size mismatch error
    tokenized_text = tokenizer(text, truncation=True, max_length=512, return_tensors="pt")

    # Perform NER using the BioBERT pipeline (truncate text)
    ner_results = ner_pipeline(tokenizer.decode(tokenized_text['input_ids'][0]), truncation=True)

    return ner_results

# Function to rejoin sub-tokens
def rejoin_subtokens(ner_results):
    combined_results = []
    current_word = ""
    current_label = None

    for entity in ner_results:
        word = entity['word'].replace("##", "")  # Rejoin sub-tokens
        label = entity['entity']

        if current_label == label:
            current_word += word
        else:
            if current_word:
                combined_results.append((current_word, current_label))
            current_word = word
            current_label = label

    if current_word:
        combined_results.append((current_word, current_label))

    return combined_results

# Function to map LABEL_0 and LABEL_1 to more meaningful labels
def map_labels(label):
    if label == "LABEL_0":
        return "Non-Entity"
    elif label == "LABEL_1":
        return "Entity"
    return label

# Function to fetch PubMed articles based on a search term
def fetch_pubmed_articles(query, num_articles, email):
    Entrez.email = email
    handle = Entrez.esearch(db="pubmed", term=query, retmax=num_articles)
    record = Entrez.read(handle)
    handle.close()

    ids = record['IdList']
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
    articles = list(Medline.parse(handle))
    handle.close()

    return articles

# Function to generate publication year distribution graph
def plot_publication_years(articles):
    years = [article.get('DP', 'No Year').split()[0] for article in articles]
    year_counts = Counter(years)

    plt.figure(figsize=(10, 6))
    plt.bar(year_counts.keys(), year_counts.values(), color='skyblue')
    plt.xlabel('Publication Year')
    plt.ylabel('Number of Articles')
    plt.title('Distribution of Articles by Publication Year')
    st.pyplot(plt)

# Function to save results to Excel
def save_to_excel(data):
    output = BytesIO()
    df = pd.DataFrame(data)
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False)
    output.seek(0)
    return output

# Streamlit UI for user inputs
st.title("PubMed Research Navigator with BioBERT NER and Data Visualization")

# Input fields for PubMed search
email = st.text_input("Enter your email (for PubMed access):")
search_term = st.text_input("Enter the PubMed search term:")
num_articles = st.number_input("Enter the number of articles to fetch:", min_value=1, max_value=1000, value=10)

# MeSH term (optional)
mesh_term = st.text_input("Enter a MeSH term (optional):")

# Fetch and analyze PubMed articles
if st.button("Fetch and Analyze PubMed Articles"):
    if email and search_term:
        st.write(f"Fetching {num_articles} articles related to '{search_term}' from PubMed...")

        # Construct PubMed query
        query = f"{search_term} AND {mesh_term}" if mesh_term else search_term
        
        # Fetch PubMed articles
        articles = fetch_pubmed_articles(query, num_articles, email)
        st.write(f"Fetched {len(articles)} articles.")

        # Extract publication year and plot graph
        st.write("Publication Year Distribution:")
        plot_publication_years(articles)

        # Process each article for NER and display results
        all_ner_results = []
        pubmed_data = []
        for article in articles:
            title = article.get('TI', 'No title available')
            abstract = article.get('AB', '')

            # Skip processing if both title and abstract are empty
            if not title.strip() and not abstract.strip():
                continue

            # Combine title and abstract for NER
            text = f"{title} {abstract}".strip()  # Ensure input text is not empty
            ner_results = bio_ner(text)
            rejoined_results = rejoin_subtokens(ner_results)

            # Store NER results with article metadata
            article_data = {
                'Title': title,
                'Abstract': abstract,
                'Publication Year': article.get('DP', 'No publication date available').split()[0],
                'Journal': article.get('TA', 'No journal available'),
                'NER Results': ", ".join([f"{word} ({map_labels(label)})" for word, label in rejoined_results])
            }
            pubmed_data.append(article_data)

            # Display NER results
            st.write(f"Title: {title}")
            for word, label in rejoined_results:
                st.write(f"Entity: {word}, Label: {map_labels(label)}")

        # Offer results as a downloadable Excel file
        excel_data = save_to_excel(pubmed_data)
        st.download_button(label="Download Excel", data=excel_data, file_name="pubmed_ner_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

    else:
        st.write("Please enter both your email and a search term.")

# Copyright information
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
