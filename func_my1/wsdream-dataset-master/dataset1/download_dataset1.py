#! /usr/bin/env python
#
# Copyright (C) 2016, WS-DREAM, CUHK
# License: MIT

import urllib
import urllib.request
import re, os
import zipfile
import shutil


url = 'https://wsdream.github.io/dataset/wsdream_dataset1'
page = urllib.request.urlopen(url)
html = page.read()
pattern = r'<a id="downloadlink" href="(.*?)"'
downloadlink = re.findall(pattern, html.decode('utf-8'))[0]
file_name = downloadlink.split('/')[-1]
page = urllib.request.urlopen(downloadlink)
meta = page.info()
file_size = int(meta.getheaders("Content-Length")[0])
print("Downloading: %s (%s bytes)" % (file_name, file_size))
urllib.urlretrieve(downloadlink, file_name)

print('Unzip data files...')
with zipfile.ZipFile(file_name, 'r') as z:
   for name in z.namelist():
      filename = os.path.basename(name)
      # skip directories
      if not filename:
         continue
      # copy file (taken from zipfile's extract)
      print(filename)
      source = z.open(name)
      target = file(filename, "wb")
      with source, target:
         shutil.copyfileobj(source, target)

os.remove(file_name)

print('==============================================')
print('Downloading data done!\n')
