#install gen3-client (installed in bin/gen3 on cloud)
wget https://github.com/uc-cdis/cdis-data-client/releases/download/2022.10/dataclient_linux.zip
unzip dataclient_linux.zip
#add to path
echo 'export PATH=$PATH:~/bin/gen3' >> ~/.bash_profile.
source ~/.bash_profile.

#configure
../bin/gen3/gen3-client configure --profile=anvil --cred=credentials.json --apiendpoint=https://gen3.theanvil.io

#download 
../bin/gen3/gen3-client download-multiple --profile=anvil --manifest=manifest.json --download-path=Crams --protocol=s3
