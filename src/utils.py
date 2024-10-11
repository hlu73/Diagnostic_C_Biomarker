# Import necessary libraries
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, RandomizedSearchCV, cross_val_score
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score, precision_recall_curve, accuracy_score, precision_score, recall_score, f1_score
from sklearn.preprocessing import StandardScaler, label_binarize, LabelEncoder
from sklearn.inspection import permutation_importance
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import learning_curve

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
    X = expression_data[features['x']]
    y = labels['Condition']
    
    return X, y

def preprocess_data(X, y):
    """
    Preprocess the data: encode target variable and split into train/test sets.
    
    Args:
    X (pd.DataFrame): Feature matrix
    y (pd.Series): Target variable
    
    Returns:
    X_train, X_test, y_train, y_test: Split datasets
    """
    # Encode the target variable
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    
    # # Split the data
    # X_train, X_test, y_train, y_test = train_test_split(X, y_encoded, test_size=0.2, random_state=42)
    
    # Use stratified sampling for train/test split
    X_train, X_test, y_train, y_test = train_test_split(X, y_encoded, test_size=0.2, random_state=42, stratify=y_encoded)
    
    return X_train, X_test, y_train, y_test

def custom_scorer(estimator, X, y):
    """
    Custom scoring function that balances accuracy and AUC-ROC.
    """
    y_pred = estimator.predict(X)
    y_pred_proba = estimator.predict_proba(X)
    
    accuracy = accuracy_score(y, y_pred)
    n_classes = len(np.unique(y))
    if n_classes > 2:
        y_bin = label_binarize(y, classes=range(n_classes))
        auc_roc = roc_auc_score(y_bin, y_pred_proba, multi_class='ovr', average='macro')
    else:
        auc_roc = roc_auc_score(y, y_pred_proba[:, 1])
    
    return (accuracy + auc_roc) / 2

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
    y_pred_proba = model.predict_proba(X_test)
    
    # Calculate evaluation metrics
    accuracy = accuracy_score(y_test, y_pred)
    
    # Handle multi-class scenario
    n_classes = len(np.unique(y_test))
    if n_classes > 2:
        # Multi-class metrics
        precision = precision_score(y_test, y_pred, average='weighted')
        recall = recall_score(y_test, y_pred, average='weighted')
        f1 = f1_score(y_test, y_pred, average='weighted')
        
        # ROC AUC for multi-class
        y_test_bin = label_binarize(y_test, classes=range(n_classes))
        auc_roc = roc_auc_score(y_test_bin, y_pred_proba, multi_class='ovr', average='macro')
    else:
        # Binary classification metrics
        precision = precision_score(y_test, y_pred)
        recall = recall_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred)
        auc_roc = roc_auc_score(y_test, y_pred_proba[:, 1])
    
    # Perform cross-validation
    cv_scores = cross_val_score(model, X_train, y_train, cv=5)
    
    evaluation_metrics = {
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1_score': f1,
        'auc_roc': auc_roc,
        'cv_scores': cv_scores,
        'cv_mean': np.mean(cv_scores),
        'cv_std': np.std(cv_scores),
        'classification_report': classification_report(y_test, y_pred)
    }
    
    feature_importance = permutation_importance(model, X_test, y_test, n_repeats=10, random_state=42)
    evaluation_metrics['feature_importance'] = feature_importance
    
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
    
# Add learning curves to detect overfitting
def plot_learning_curve(estimator, X, y, cv=None, n_jobs=None, train_sizes=np.linspace(.1, 1.0, 5)):
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=cv, n_jobs=n_jobs, train_sizes=train_sizes, scoring=custom_scorer)
    
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    
    plt.figure()
    plt.title("Learning Curve")
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    plt.grid()
    
    plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1, color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r", label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g", label="Cross-validation score")
    
    plt.legend(loc="best")
    plt.show()