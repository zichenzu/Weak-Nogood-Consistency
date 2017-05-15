# Weak-Nogood-Consistency
Weak Nogood Consistency

This is the implementation of the method LReSBDS using the incNGs global constraint proposed in the paper "Jimmy H.M. Lee and Zichen Zhu.  Filtering Nogoods Lazily in Dynamic Symmetry Breaking During Search, Proceedings of the 24th International Joint Conference on Artificial Intelligence (IJCAI 2015), pages 339-345, Buenos Aires, Argentina, July, 2015"

Please use Gecode Solver 4.2.0 to run these files.

Put WNC folder into gecode folder and efpa_lresbds_WNC.cpp file into example folder.

To run the N-Queens problem (e.g. N=10) using one of the symmetry breaking methods NONE/DoubleLex/SnakeLex/LDSB/LReSBDS:

./efpa_lresbds_WNC -search none 5 3 3 4

./efpa_lresbds_WNC -search double 5 3 3 4

./efpa_lresbds_WNC -search snake 5 3 3 4

./efpa_lresbds_WNC -search ldsb 5 3 3 4

./efpa_lresbds_WNC -search lresbds 5 3 3 4
