# Combine-minimizer
The  source code of combine minimizer for sequence-to-graph alignment.

Combine minimizer is the combined hashing of adjacent minimizer seeds. The main work  borrowed giraffe's alignment framework and modified the seed selection and positioning parts(vg/deps/include/gbwtgraph).
If run the code, change the include files.

And change find-minimizers function in minimizer_mapper.cpp
```
 vector<tuple<gbwtgraph::DefaultMinimizerIndex::minimizer_type, size_t, size_t>> minimizers =
            this->minimizer_index.minimizer_regions_combine_gap(sequence, this->minimizer_index.w());
```


