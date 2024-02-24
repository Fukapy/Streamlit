import streamlit as st
from chempy import Substance

def calculate_molecular_weight(formula):
    try:
        substance = Substance.from_formula(formula)
        molecular_weight = substance.mass
        return molecular_weight
    except Exception as e:
        return str(e)

def main():
    st.title("分子式から分子量を算出")

    # ユーザーに分子式を入力させる
    formula_input = st.text_input("分子式を入力してください (例: H2O, C6H12O6)")

    # 分子量の計算
    if st.button("計算"):
        if formula_input:
            molecular_weight = calculate_molecular_weight(formula_input)
            st.success(f"{formula_input} の分子量は {molecular_weight:.2f} g/mol です。")
        else:
            st.warning("分子式を入力してください。")

if __name__ == "__main__":
    main()
