# Merge 以下の説明

- DB ... **D**ESHIMA **D**ata**B**ase (DDB) 関連
  - DDB_20171018151938.fits
    - DDB (FITS 形式) KIDSINFO と RCP 情報が入っている
  - DDB_dict.yaml
    - DDB------.fits を作る際の各HDU の Key や Comment 等が保存されているファイル
  - create_DDB.ipynb
    - DDB------.fits を作るコード。``DDB_dict.yaml`` を用いる
  - create_DDBDict.ipynb
    - ``DDB_dict.yaml`` を作るためのコード (jupyter notebook を用いる)


- datafile ... DFITS を作る際に使用するデータファイル (antennalog etc.)


- DFITS_dict.yaml
  - DFITS-----.fits を作る際の各HDU の Key や Comment 等が保存されているファイル
- DFITS_merge.py
  - DFITS-----.fits のための HDUList を返すスクリプト
- create_DFITSDict.ipynb
  - ``DFITS_dict.yaml`` を作るコード
- sample_dfits_fromaste.ipynb
  - DFITS_merge.py の使用例
