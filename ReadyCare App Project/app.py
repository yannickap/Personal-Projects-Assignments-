import streamlit as st
import pandas as pd

st.set_page_config(page_title="ReadyCare Telehealth Readiness", layout="wide")

st.title("📱 ReadyCare: Telehealth Readiness Tool")

st.write("Answer the questions below to check your digital readiness for telehealth visits.")

questions = {
    "q1": "Do you feel confident using a smartphone or tablet?",
    "q2": "Can you join a video call without help?",
    "q3": "Do you know how to check your internet connection?",
    "q4": "Do you feel comfortable downloading/opening apps?",
    "q5": "Do you know how to adjust microphone settings?",
    "q6": "Do you know how to turn your camera on/off?",
    "q7": "Have you attended a telehealth appointment before?",
    "q8": "Do you have a private, quiet space?",
    "q9": "Do you have a working camera?",
    "q10": "Do you have reliable internet?"
}

options = {
    "Yes, without help": 2,
    "Yes, but sometimes need help": 1,
    "No": 0
}

scores = []

for key, q in questions.items():
    choice = st.radio(q, list(options.keys()), key=key)
    scores.append(options[choice])

if st.button("Calculate My Readiness"):
    total = sum(scores)
    max_score = len(scores) * 2
    pct = (total / max_score) * 100

    if pct >= 75:
        level = "HIGH"
        color = "green"
    elif pct >= 40:
        level = "MEDIUM"
        color = "orange"
    else:
        level = "LOW"
        color = "red"
    
    st.markdown(f"## **Your Score: {pct:.1f}% — <span style='color:{color}'>{level} readiness</span>**", unsafe_allow_html=True)

    st.subheader("Recommendations")
    
    if pct < 75:
        st.write("- Practice joining video calls using Zoom or FaceTime.")
    if pct < 60:
        st.write("- Ask your clinic about low-cost internet programs.")
    if pct < 40:
        st.write("- Request a telehealth navigator from your clinic.")
