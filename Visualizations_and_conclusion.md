### Visualizations

Here are some key plots illustrating the findings:

#### 1. Expression of Top Stable Genes by Disease Status

This visualization is critical as it reveals the distinct expression patterns of the genes identified through stability selection. Notice the clear separation in expression levels for each gene between the 'Disease' and 'Normal' groups, which explains the observed perfect classification performance.

<img width="913" height="896" alt="Expression_of _top_stable_genes" src="https://github.com/user-attachments/assets/a8103a02-569e-4c82-9bd3-6be81dae6d4b" />

*(This plot highlights how these selected genes act as perfect biomarkers within this dataset.)*

#### 3. Random Forest ROC Curves (Training vs. Test Set)

This plot shows the Receiver Operating Characteristic (ROC) curves for the Random Forest model on both the training and independent test sets. An AUC of 1.0 signifies perfect classification, where the model achieves a 100% true positive rate without any false positives. The overlap of the training and test curves further supports that this perfect performance generalizes to unseen data *within this dataset*.

<img width="821" height="896" alt="c30d72f0-ce87-4768-b14a-cc9e84d64068" src="https://github.com/user-attachments/assets/2bd906bc-c007-4caa-a79d-f8b7a7827c44" />

*(The LASSO ROC curve would look identical in this scenario, also indicating an AUC of 1.0.)*

---

## Conclusion & Future Directions

This robust machine learning pipeline has identified a strong gene signature capable of perfectly distinguishing IPF patients from healthy controls within this dataset. The consistent AUC of 1.0, even on an independent test set, is directly attributable to the absolute separation in expression levels of the identified stable genes.
