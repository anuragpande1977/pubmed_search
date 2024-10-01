import os
import pandas as pd
import streamlit as st
from Bio import Entrez, Medline
from io import BytesIO
import matplotlib.pyplot as plt

# Define PubMed article types and their corresponding search tags
article_types = {
    "Clinical Trials": "Clinical Trial[pt]",
    "Meta-Analysis": "Meta-Analysis[pt]",
    "Randomized Controlled Trials": "Randomized Controlled Trial[pt]",
    "Reviews": "Review[pt]",
    "Systematic Reviews": "Systematic Review[pt]",
}

# Function to construct query
def construct_query(search_term, choice):
    chosen_article_type = article_types[choice]
    query = f"({search_term}) AND {chosen_article_type}"
    return query

# Function to fetch articles from PubMed
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

# Function to generate Excel file in memory
def save_to_excel(articles):
    output = BytesIO()
    
    data = [{
        'Title': article.get('TI', 'No title available'),
        'Authors': ', '.join(article.get('AU', 'No authors available')),
        'Abstract': article.get('AB', 'No abstract available'),
        'Publication Date': article.get('DP', 'No publication date available'),
        'Journal': article.get('TA', 'No journal available')
    } for article in articles]
    
    df = pd.DataFrame(data)
    
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False)
    
    output.seek(0)  # Move the buffer position to the beginning
    return output

# Function to count articles by publication year
def count_articles_by_year(articles):
    year_count = {}
    
    for article in articles:
        pub_date = article.get('DP', 'No publication date available')
        if pub_date != 'No publication date available':
            year = pub_date.split()[0]  # Extracting the year
            if year.isdigit():  # Ensure it's a valid year
                year_count[year] = year_count.get(year, 0) + 1
    
    return year_count

# Function to display the bar chart
def plot_publication_years(year_count):
    if year_count:
        years = list(year_count.keys())
        counts = list(year_count.values())

        plt.figure(figsize=(10, 5))
        plt.bar(years, counts, color='blue')
        plt.xlabel('Publication Year')
        plt.ylabel('Number of Articles')
        plt.title('Number of Articles Published per Year')
        plt.xticks(rotation=45)
        st.pyplot(plt)

# Streamlit UI for user inputs
st.title("PubMed Article Fetcher")
st.write("Search PubMed for articles and save the results as an Excel file.")

email = st.text_input("Enter your email (for PubMed access):")
search_term = st.text_input("Enter the general search term:")
article_choice = st.selectbox("Select article type:", list(article_types.keys()))
num_articles = st.number_input("Enter the number of articles to fetch:", min_value=1, max_value=1000, value=10)

if st.button("Fetch Articles"):
    if email and search_term:
        query = construct_query(search_term, article_choice)
        articles = fetch_abstracts(query, num_articles, email)
        
        if articles:
            excel_data = save_to_excel(articles)
            st.download_button(
                label="Download Excel file",
                data=excel_data,
                file_name="pubmed_articles.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
            
            # Count articles by year and plot the results
            year_count = count_articles_by_year(articles)
            if year_count:
                st.write("Publication Year Distribution:")
                plot_publication_years(year_count)
            else:
                st.write("No valid publication dates found to plot the graph.")
        else:
            st.write("No articles fetched.")
    else:
        st.write("Please fill in all the required fields.")

# Copyright information
st.text("Software copyright holder: Anurag Pande")
