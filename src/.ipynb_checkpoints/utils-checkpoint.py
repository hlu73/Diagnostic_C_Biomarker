# Import necessary libraries
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, RandomizedSearchCV, cross_val_score
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_auc_score, precision_recall_curve
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# Utility functions

def load_data(expression_path, labels_path, feature_set_path):
    """
    Load expression data, labels, and feature set from CSV files.
    
    Args:
    expression_path (str): Path to the CSV file containing gene expression data
    labels_path (str): Path to the CSV file containing labels
    feature_set_path (str): Path to the CSV file containing selected features
    
    Returns:
    X (pd.DataFrame): Feature matrix
    y (pd.Series): Target variable
    """
    # Load the expression data
    expression_data = pd.read_csv(expression_path)
    
    # Load the labels
    labels = pd.read_csv(labels_path)
    
    # Load the selected features
    features = pd.read_csv(feature_set_path)
    
    # Extract features and target variable
    X = expression_data[features['gene_name']]
    y = labels['Condition']
    
    return X, y

def preprocess_data(X, y):
    """
    Preprocess the data: split into train/test sets and scale features.
    
    Args:
    X (pd.DataFrame): Feature matrix
    y (pd.Series): Target variable
    
    Returns:
    X_train, X_test, y_train, y_test: Split and scaled datasets
    """
    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Scale the features
    #scaler = StandardScaler()
    #X_train_scaled = scaler.fit_transform(X_train)
    #X_test_scaled = scaler.transform(X_test)
    
    return X_train, X_test, y_train, y_test

def train_and_evaluate_model(model, X_train, X_test, y_train, y_test):
    """
    Train a model and evaluate its performance.
    
    Args:
    model: The machine learning model to train
    X_train, X_test, y_train, y_test: Training and testing data
    
    Returns:
    trained_model: The trained model
    evaluation_metrics: Dictionary containing various evaluation metrics
    """
    # Train the model
    model.fit(X_train, y_train)
    
    # Make predictions
    y_pred = model.predict(X_test)
    y_pred_proba = model.predict_proba(X_test)[:, 1]
    
    # Calculate evaluation metrics
    accuracy = accuracy_score(y_test, y_pred)
    auc_roc = roc_auc_score(y_test, y_pred_proba)
    precision, recall, _ = precision_recall_curve(y_test, y_pred_proba)
    
    # Perform cross-validation
    cv_scores = cross_val_score(model, X_train, y_train, cv=5)
    
    evaluation_metrics = {
        'accuracy': accuracy,
        'auc_roc': auc_roc,
        'precision': precision,
        'recall': recall,
        'cv_scores': cv_scores,
        'classification_report': classification_report(y_test, y_pred)
    }
    
    return model, evaluation_metrics

def plot_feature_importance(model, feature_names):
    """
    Plot feature importance for a given model.
    
    Args:
    model: Trained model with feature_importances_ attribute
    feature_names: List of feature names
    """
    feature_importance = pd.DataFrame({
        'feature': feature_names,
        'importance': model.feature_importances_
    }).sort_values('importance', ascending=False)

    plt.figure(figsize=(10, 6))
    sns.barplot(x='importance', y='feature', data=feature_importance.head(20))
    plt.title('Top 20 Most Important Features')
    plt.tight_layout()
    plt.show()

def plot_confusion_matrix(y_true, y_pred):
    """
    Plot confusion matrix.
    
    Args:
    y_true: True labels
    y_pred: Predicted labels
    """
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.title('Confusion Matrix')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    plt.show()