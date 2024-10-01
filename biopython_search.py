import pandas as pd
import streamlit as st
from transformers import AutoTokenizer, AutoModelForTokenClassification, pipeline

# Load BioBERT Tokenizer and Model for NER
tokenizer = AutoTokenizer.from_pretrained("dmis-lab/biobert-base-cased-v1.1")
model = AutoModelForTokenClassification.from_pretrained("dmis-lab/biobert-base-cased-v1.1")

# Initialize NER pipeline using BioBERT
ner_pipeline = pipeline("ner", model=model, tokenizer=tokenizer)

# Function to perform NER on biomedical text using BioBERT
def bio_ner(text):
    ner_results = ner_pipeline(text)
    return ner_results

# Streamlit UI for user inputs
st.title("PubMed Research Navigator with BioBERT NER")
email = st.text_input("Enter your email (for PubMed access):")
search_term = st.text_input("Enter the general search term:")
num_articles = st.number_input("Enter the number of articles to fetch:", min_value=1, max_value=1000, value=10)

# Sample biomedical text (this would come from fetched PubMed articles)
sample_text = "SARS-CoV-2, the virus responsible for COVID-19, affects the respiratory system."

if st.button("Analyze Text with BioBERT NER"):
    if email and search_term:
        st.write("Performing NER with BioBERT...")
        
        # Call the BioBERT NER function on the sample text
        ner_results = bio_ner(sample_text)
        
        # Display results
        if ner_results:
            st.write("Named Entities extracted from the text:")
            for entity in ner_results:
                st.write(f"Entity: {entity['word']}, Label: {entity['entity']}, Score: {entity['score']}")
        else:
            st.write("No entities found.")
    else:
        st.write("Please enter valid input.")

# Display copyright information
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
