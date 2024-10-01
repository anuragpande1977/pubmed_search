import os
import pandas as pd
import streamlit as st
from Bio import Entrez, Medline
from io import BytesIO
import matplotlib.pyplot as plt
from pyvis.network import Network
import ast

# Define PubMed article types and their corresponding search tags
article_types = {
    "Clinical Trials": "Clinical Trial[pt]",
    "Meta-Analysis": "Meta-Analysis[pt]",
    "Randomized Controlled Trials": "Randomized Controlled Trial[pt]",
    "Reviews": "Review[pt]",
    "Systematic Reviews": "Systematic Review[pt]",
}

# Function to construct query with optional MeSH term
def construct_query(search_term, mesh_term, choice):
    chosen_article_type = article_types[choice]
    query = f"({search_term}) AND {chosen_article_type}"
    
    # Include MeSH term if provided
    if mesh_term:
        query += f" AND {mesh_term}[MeSH Terms]"
    
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

# Function to extract MeSH terms and Chemicals from articles
def extract_entities(articles):
    extracted_data = []
    for article in articles:
        title = article.get('TI', 'No title available')
        mesh_terms = article.get('MH', [])
        abstract = article.get('AB', 'No abstract available')

        chemicals = []  # You might extract chemicals from 'RN' or 'OT' fields or using custom extraction logic
        # Simplified chemical extraction logic - you might want to use NLP or regex in a real application
        chemicals += [word for word in abstract.split() if 'chem' in word.lower()]  # Example chemical detection

        extracted_data.append({
            'Title': title,
            'MeSH Terms': mesh_terms,
            'Chemicals': chemicals
        })
    
    return pd.DataFrame(extracted_data)

# Function to create interactive graph between MeSH terms and Chemicals
def create_interactive_graph(df, title_column, mesh_column, chemical_column, central_word, background_color, output_file_path):
    color_options = {
        "mesh": "#FF5733",  # Example color for MeSH terms
        "chemical": "#33CFFF",  # Example color for chemicals
    }
    default_node_size = 15  # Uniform node size for all nodes
    central_color = "#FFFF00"  # Highlight color for the central entity

    net = Network(notebook=False, height="100%", width="100%", bgcolor=background_color, font_color="white")
    net.barnes_hut()

    # Explicitly add the central entity
    net.add_node(central_word, label=central_word, color=central_color, size=default_node_size, title="Central Entity")

    for _, row in df.iterrows():
        mesh_terms = row[mesh_column]
        chemicals = row[chemical_column]
        title = row[title_column]

        for term in mesh_terms:
            if not net.get_node(term):
                net.add_node(term, label=term, color=color_options["mesh"], size=default_node_size, title=title)
            net.add_edge(central_word, term)

        for chemical in chemicals:
            if not net.get_node(chemical):
                net.add_node(chemical, label=chemical, color=color_options["chemical"], size=default_node_size, title=title)
            net.add_edge(central_word, chemical)

    net.write_html(output_file_path)
    return output_file_path

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

# Streamlit UI for user inputs
st.title("PubMed Research Navigator with MeSH-Chemical Relationship Graph")
st.write("Search PubMed for articles, save results as an Excel file, or create a relationship graph between MeSH terms and chemicals.")

email = st.text_input("Enter your email (for PubMed access):")
search_term = st.text_input("Enter the general search term:")
mesh_term = st.text_input("Enter an optional MeSH term (leave blank if not needed):")
article_choice = st.selectbox("Select article type:", list(article_types.keys()))
num_articles = st.number_input("Enter the number of articles to fetch:", min_value=1, max_value=1000, value=10)
central_word = st.text_input("Enter the central word for the graph (e.g., disease or keyword):").lower()

if st.button("Fetch Articles and Create Graph"):
    if email and search_term and central_word:
        query = construct_query(search_term, mesh_term, article_choice)
        articles = fetch_abstracts(query, num_articles, email)
        
        if articles:
            # Save articles as an Excel file
            excel_data = save_to_excel(articles)
            st.download_button(
                label="Download Excel file",
                data=excel_data,
                file_name="pubmed_articles.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

            # Extract entities (MeSH terms and chemicals)
            df_entities = extract_entities(articles)

            # Create and visualize the interactive graph
            background_color = 'black'
            output_html = 'mesh_chemical_relationship_graph.html'
            graph_path = create_interactive_graph(df_entities, 'Title', 'MeSH Terms', 'Chemicals', central_word, background_color, output_html)
            st.write(f"Graph created! You can download or view it [here](./{output_html}).")
        else:
            st.write("No articles fetched.")
    else:
        st.write("Please fill in all the required fields.")

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
