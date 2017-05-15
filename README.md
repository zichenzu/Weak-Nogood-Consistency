# Weak-Nogood-Consistency
Weak Nogood Consistency

This is the implementation of the Weak-Nogood-Consistency (WNC) on the symmetry breaking nogoods proposed in the paper "Jimmy H.M. Lee and Zichen Zhu.  Filtering Nogoods Lazily in Dynamic Symmetry Breaking During Search, Proceedings of the 24th International Joint Conference on Artificial Intelligence (IJCAI 2015), pages 339-345, Buenos Aires, Argentina, July, 2015"

Please use Gecode Solver 4.2.0 to run these files.

Put WNC folder into gecode folder and efpa_lresbds_WNC.cpp file into example folder.

To run the EFPA problem (e.g. (3 4 6 4)) using one of the symmetry breaking methods NONE/DoubleLex/SnakeLex/LDSB/LReSBDS:

./efpa_lresbds_WNC -search none 3 4 6 4

./efpa_lresbds_WNC -search double 3 4 6 4

./efpa_lresbds_WNC -search snake 3 4 6 4

./efpa_lresbds_WNC -search ldsb 3 4 6 4

./efpa_lresbds_WNC -search lresbds 3 4 6 4
