# GROMACS 編譯方法

###### tags: `Linux` `計算化學` `NCHC`

## 前言

GROMACS為開源的化學計算軟體，以下文章為NCHC安裝軟體之用，內文皆參考自下列官方網址:
- [GROMACS官方網站](https://manual.gromacs.org/2024.1/index.html)
- [Github](https://github.com/gromacs/gromacs)
:::warning
請注意，目前文章為事前準備之用，無實際上機操作，因此內容可能有錯誤之虞。
:::

## 環境配置

目前計畫先行於T3上進行安裝，故以下安裝環境以T3之**Intel**及**AMD**節點為考量，並依節點之處理器不同採不同之配置以確保最優化之效能。同時T3上**無配置GPU**，因此此文章編譯時就先不考慮軟體的GPU支援了。

- **Intel節點**

| 項目 | 配置 |
| :--: | :--: |
| Compiler | oneAPI |
| MPI | OpenMPI |

- **AMD節點**

| 項目 | 配置 |
| :--: | :--: |
| Compiler | GUN |
| MPI | OpenMPI |

## [編譯方法](https://manual.gromacs.org/2024.1/install-guide/index.html#cmake-advanced-options)

1. **下載並解壓縮檔案後進入目錄**

```bash!
wget https://ftp.gromacs.org/gromacs/gromacs-2024.1.tar.gz && tar xfz gromacs-2024.1.tar.gz && cd gromacs* && ls && mkdir build && cd build
```

2. **載入編譯用模組**

- :::spoiler Intel節點
  ```bash!
    ml purge && ml load compiler/intel/2020u4 openmpi/4.1.1 tools/cmake/3.21.2
  ```
  :::
  
- :::spoiler AMD節點
  ```bash!
    ml purge && ml load compiler/gcc/11.2.0 OpenMPI/4.1.4 tools/cmake/3.21.2
  ```
  :::
  
3. **使用cmake產生建構檔**

cmake會自動使用環境中載入的編譯器進行編譯，因此安裝時不用特別指定編譯器。

```bash!
cmake .. -DGMX_MPI=on -DCMAKE_PREFIX_PATH=</path/to/needs/libraries> -DCMAKE_INSTALL_PREFIX=</path/where/gromacs/is/installed> -DGMX_BUILD_OWN_FFTW=ON
```

以下依添加的option逐個說明:
- **-DGMX_MPI=** 開啟MPI支援
- **-DCMAKE_PREFIX_PATH=** 增加庫的路徑
- **-DCMAKE_INSTALL_PREFIX=** 自定義GROMACS的安裝路徑
- **-DGMX_BUILD_OWN_FFTW=** 允許GROMACS建立自己的FFTW

4. **編譯並安裝GROMACS**

```bash!=
make
make check
sudo make install
source /usr/local/gromacs/bin/GMXRC
```

## 測試(以Intel節點為例)

[測試檔案及教學來源](http://www.mdtutorials.com/gmx/lysozyme/index.html)，對指令內容有相關問題請參考連結內說明。

:::spoiler {state="open"} **水中溶菌酶模擬**

1. 下載蛋白質pdb檔並去除水分子

```bash!
wget https://files.rcsb.org/download/1AKI.pdb && grep -v HOH 1aki.pdb > 1AKI_clean.pdb
```

2. 產生拓譜結構、位置結構等需要的結構檔

```bash!
gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce
```

3. 選擇力場

```bash!
>15
```

4. 定義cell大小

```bash!
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic
```

5. 在cell中填滿溶液

```bash!
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top
```

6. 下載後續處理需要的動力學參數檔

```bash!
wget http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp
```

7. 將動力學參數檔與拓譜資訊結合產生出運行輸入文件(替換水分子以達系統電中性)

```bash!
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr
```

8. 將文件傳遞給genion

```bash!
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```

9. 將檔案組合轉換成二進制檔

```bash!
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
```

10. 執行能量最小化(結構優化)

```bash!
gmx mdrun -v -deffnm em
```

11. 分析執行最小化後的能量曲線，確認是否收斂

```bash!
gmx energy -f em.edr -o potential.xvg
```

12. 下載參數檔並執行等溫等容的平衡計算

```bash!
wget http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt
```

13. 分析計算結果的溫度變化

```bash!
gmx energy -f nvt.edr -o temperature.xvg

>16 0
```

14. 下載參數檔並執行等溫等壓的平衡計算

```bash!
wget http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun-deffnm npt
```

15. 分析計算結果的壓力變化

```bash!
gmx energy -f npt.edr -o pressure.xvg

>18 0
```
15. 分析計算結果的密度變化

```bash!
gmx energy -f npt.edr -o density.xvg

>24 0
```

16. 下載腳本並放開系統開始進行動力學模擬計算

```bash!
wget http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1
```

17. 使用trjconv工具處理計算後的結果供後續分析使用

```bash!
gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
```

18. 分析RMSD

```bash!
gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
```





:::


:::spoiler {state="open"} **JOB腳本與執行範例**

```bash!
#!/bin/bash
#SBATCH -A XXXXXXXXX        # Project ID
#SBATCH -J test             # Job name
#SBATCH -p test             # Partition name
#SBATCH -n 24               # Number of MPI tasks (i.e. processes)
#SBATCH -c 1                # Number of cores per MPI task
#SBATCH -N 3                # Maximum number of nodes to be allocated
#SBATCH -o %j.out           # Path to the standard output file
#SBATCH -e %j.err           # Path to the standard error ouput file

module purge
ml compiler/intel/2020u4 openmpi/4.1.1

mpirun /The/pathway/you/install/GROMACS ...
```

```bash!
sbatch test.sh
```

:::

