import os
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from io import BytesIO
from Bio import Entrez, Medline
from collections import Counter

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

# Function to generate a pie chart for publication year distribution with number of articles
def plot_publication_years_pie(articles):
    years = [article.get('DP', 'No Year').split()[0] for article in articles]
    year_counts = Counter(years)

    plt.figure(figsize=(7, 7))
    colors = plt.cm.Set3.colors  # Use a bright color palette
    plt.pie(
        year_counts.values(), 
        labels=[f"{year} ({count})" for year, count in year_counts.items()],  # Show the number of articles
        startangle=140, 
        colors=colors
    )
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
st.title("PubMed Research Navigator with Data Visualization")

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

        # Prepare PubMed data for download
        pubmed_data = []
        for article in articles:
            article_data = {
                'Title': article.get('TI', 'No title available'),
                'Abstract': article.get('AB', 'No abstract available'),
                'Publication Year': article.get('DP', 'No publication date available').split()[0],
                'Journal': article.get('TA', 'No journal available')
            }
            pubmed_data.append(article_data)

        # Offer results as a downloadable Excel file
        excel_data = save_to_excel(pubmed_data)
        st.download_button(label="Download Excel", data=excel_data, file_name="pubmed_articles.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

        # Extract publication year and plot pie chart
        st.write("Publication Year Distribution:")
        plot_publication_years_pie(articles)

    else:
        st.write("Please enter both your email and a search term.")

# Copyright information
st.write("""
    ### Copyright Information
    Copyright (c) 2024 Anurag Pande
    
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
to use the Software solely for non-commercial and research purposes, provided that the following conditions are met:

1. The Software and its associated documentation must not be modified or redistributed in any form.
2. The Software may not be used for commercial purposes or for any purpose that generates profit.
3. Any use of the Software must include the original copyright notice and this permission notice in all copies or substantial portions of the Software.
4. The results produced by the Software are provided for research purposes only and should not be relied upon for any clinical or legal decision-making.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
IN THE SOFTWARE.

By using this Software, you agree to the above terms and conditions.

""")
