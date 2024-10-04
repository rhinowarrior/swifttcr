for file in *_merged.pdb; do
    grep '^ATOM' "$file" > temp && mv temp "$file"
done