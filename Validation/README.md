
python Model/TrunKitten/pipeline/predict.py \
  --annotated output/LMNA_TrunKitten_10features_with_id.tsv \
  --model Model/TrunKitten/model/trunkitten.pkl \
  --metadata Model/TrunKitten/model/trunkitten_features.json \
  --out output/LMNA_TrunKitten_predictions.tsv
  
