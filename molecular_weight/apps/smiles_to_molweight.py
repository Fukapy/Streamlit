import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_molecular_weight(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        molecular_weight = Descriptors.MolWt(mol)
        return molecular_weight
    except Exception as e:
        return str(e)

def main():
    st.title("SMILESから分子量を算出")

    # ユーザーにSMILES表記を入力させる
    smiles_input = st.text_input("SMILES表記を入力してください (例: c1ccccc1, C#N)")

    # 分子量の計算
    if st.button("計算"):
        if smiles_input:
            molecular_weight = calculate_molecular_weight(smiles_input)
            st.success(f"分子量は {molecular_weight:.2f} g/molです。")
        else:
            st.warning("SMILES表記を入力してください。")

if __name__ == "__main__":
    main()
