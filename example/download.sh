wget --no-check-certificate https://figshare.com/ndownloader/files/44813476 -O data.zip &&
unzip data.zip &&
mv example_data/* . &&
rm -r data.zip example_data
