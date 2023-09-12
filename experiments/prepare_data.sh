mkdir datasets
wget -c https://www.openml.org/data/download/1593761/phpTXWrKb -O datasets/heart.arff
wget -c https://www.openml.org/data/download/44/dataset_44_spambase.arff -O datasets/spam.arff
wget -c https://www.openml.org/data/download/1592296/php9xWOpn -O datasets/steel.arff
wget -c https://www.openml.org/data/download/1592318/phpAmSP4g -O datasets/breast.arff
mkdir datasets_quantized
mkdir datasets_synthetic
python3 generate_trees.py # Generates all necessary decision trees