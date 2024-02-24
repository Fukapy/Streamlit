import streamlit as st
from multiapp import MultiApp
from apps import molform_to_molweight, smiles_to_molweight, smiles_to_molform
app = MultiApp() 

app.add_app("分子式から分子量を算出", molform_to_molweight.main)
app.add_app("SMILESから分子量を算出", smiles_to_molweight.main) 
app.add_app("SMILESから分子式に変換", smiles_to_molform.main) 

app.run()