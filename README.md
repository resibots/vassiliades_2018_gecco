# Source code for Discovering the Elite Hypervolume by Leveraging Interspecies Correlation (GECCO 2018)

V. Vassiliades & J.-B. Mouret

## How to properly clone this repo

```
git clone --recursive https://github.com/resibots/vassiliades_2018_gecco.git
```

## Dependencies

- [Boost]: C++ Template Libraries
- [Eigen]: Linear Algebra C++ Library
- [realpath]: `sudo apt-get install realpath`

## How to easily compile everything

**Important:** Make sure you have installed all the dependencies of each repo. Otherwise the build will fail.

From the root of this repo run:

```
sh compile_all.sh
```

## How to run the experiments
```
cd code/sferes2/build/exp/gecco2018exp/src/experiments/
```

### Schwefel function
```
./schwefel_variation<XXX> ../../../../../exp/gecco2018exp/centroids/centroids_10000_2_normalized_-5.0_5.0.dat
```

### Arm  
```
./arm_variation<XXX> ../../../../../exp/gecco2018exp/centroids/centroids_10000_2_normalized_-1.0_1.0.dat
```

### Hexapod
```
./hexa_variation<XXX> ../../../../../exp/gecco2018exp/centroids/centroids_10000_6.dat
```

### Replace `<XXX>` with:

- `isolinedd` for the **Iso+LineDD** operator
- `linedd` for the **LineDD** operator
- `iso` for the **Iso** operator
- `isodd` for the **IsoDD** operator
- `isosa` for the **IsoSA** operator
- `gc` for the **Global Correlation** operator
- `sbx` for the **SBX** operator


## How to easily clean everything

From the root of this repo run:

```
sh clear_all.sh
```

## LICENSE

[CeCILL]

[CeCILL]: http://www.cecill.info/index.en.html
[Boost]: http://www.boost.org
[Eigen]: http://eigen.tuxfamily.org/
[realpath]: http://manpages.ubuntu.com/manpages/jaunty/man1/realpath.1.html
