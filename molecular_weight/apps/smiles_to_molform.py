import streamlit as st
from rdkit import Chem

def convert_smiles_to_formula(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        formula = calc_formula_with_hydrogen(mol)
        return formula
    except Exception as e:
        return str(e)

def calc_formula_with_hydrogen(mol):
    # 分子式の計算（水素も含む）
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    hydrogen_count = sum([a.GetTotalNumHs() for a in mol.GetAtoms()])  # 各原子に結合している水素の数を合計
    
    formula = ""
    
    for atom in set(atoms):
        count = atoms.count(atom)
        if count > 1:
            formula += f"{atom}{count}"
        else:
            formula += atom
    
    if hydrogen_count > 0:
        formula += f"H{hydrogen_count}"
    
    return formula

def main():
    st.title("SMILESから分子式に変換")

    # ユーザーにSMILES表記を入力させる
    smiles_input = st.text_input("SMILES表記を入力してください (例: c1ccccc1, C#N)")

    # 分子式への変換
    if st.button("変換"):
        if smiles_input:
            molecular_formula = convert_smiles_to_formula(smiles_input)
            st.success(f"SMILES表記 {smiles_input} の分子式は {molecular_formula} です。")
        else:
            st.warning("SMILES表記を入力してください。")

if __name__ == "__main__":
    main()
