import os
import pandas as pd
import streamlit as st
from Bio import Entrez, Medline

# Define PubMed article types and their corresponding search tags
article_types = {
    "Clinical Trials": "Clinical Trial[pt]",
    "Meta-Analysis": "Meta-Analysis[pt]",
    "Randomized Controlled Trials": "Randomized Controlled Trial[pt]",
    "Reviews": "Review[pt]",
    "Systematic Reviews": "Systematic Review[pt]",
}

# Streamlit UI for user inputs
st.title("PubMed Article Fetcher")
st.write("Search PubMed for articles and save the results as an Excel file.")

email = st.text_input("Enter your email (for PubMed access):")
search_term = st.text_input("Enter the general search term:")
article_choice = st.selectbox("Select article type:", list(article_types.keys()))
num_articles = st.number_input("Enter the number of articles to fetch:", min_value=1, max_value=1000, value=10)
directory = st.text_input("Enter the directory to save the file (default is current directory):", ".")
filename = st.text_input("Enter the filename for the Excel file (with .xlsx extension):", "results.xlsx")

def construct_query(search_term, choice):
    chosen_article_type = article_types[choice]
    query = f"({search_term}) AND {chosen_article_type}"
    return query

def fetch_abstracts(query, num_articles, email):
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

        return articles
    except Exception as e:
        st.write(f"An error occurred: {e}")
        return []

def save_to_excel(articles, directory, filename):
    if not articles:
        st.write("No articles to save.")
        return
    full_path = os.path.join(directory, filename)
    data = [{
        'Title': article.get('TI', 'No title available'),
        'Authors': ', '.join(article.get('AU', 'No authors available')),
        'Abstract': article.get('AB', 'No abstract available'),
        'Publication Date': article.get('DP', 'No publication date available'),
        'Journal': article.get('TA', 'No journal available')
    } for article in articles]
    
    df = pd.DataFrame(data)
    if not os.path.exists(directory):
        os.makedirs(directory)
    df.to_excel(full_path, index=False, engine='openpyxl')
    st.write(f"Data saved to {full_path}")

# When the user clicks the 'Fetch Articles' button
if st.button("Fetch Articles"):
    if email and search_term and filename:
        query = construct_query(search_term, article_choice)
        articles = fetch_abstracts(query, int(num_articles), email)
        if articles:
            save_to_excel(articles, directory, filename)
        else:
            st.write("No articles fetched.")
    else:
        st.write("Please fill in all the required fields.")
