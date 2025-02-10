
import sys
import argparse
sys.path.append(".")
sys.path.append("/home/share/huadjyin/home/zhoumin3/LucaOne/LucaOneApp-master")
sys.path.append("/home/share/huadjyin/home/zhoumin3/LucaOne/LucaOneApp-master/algorithms")
try:
    from algorithms.llm.esm.predict_embedding import main
except ImportError:
    from algorithms.llm.esm.predict_embedding import main

def get_args():
    parser = argparse.ArgumentParser(description='ESM2 Embedding')
    # for one seq
    parser.add_argument("--seq_id", type=str, default=None,
                        help="the seq id")
    parser.add_argument("--seq", type=str, default=None,
                        help="when to input a seq")
    parser.add_argument("--seq_type", type=str, default="prot",
                        choices=["prot"],
                        help="the input seq type")

    # for many
    parser.add_argument("--input_file", type=str, default=None,
                        help="the input filepath(.fasta or .csv or .tsv)")

    # for input csv/tsv
    parser.add_argument("--id_idx", type=int, default=None,
                        help="id col idx(0 start)")
    parser.add_argument("--seq_idx", type=int, default=None,
                        help="seq col idx(0 start)")

    # for saved path
    parser.add_argument("--save_path", type=str, default=None,
                        help="embedding file save dir path")

    parser.add_argument("--llm_type", type=str, default="esm2",
                        choices=["esm2", "esm", "ESM"],
                        help="llm type")
    parser.add_argument("--llm_version", type=str, default="3B",
                        choices=["15B", "3B", "650M", "150M"],
                        help="llm version")

    # for embedding
    parser.add_argument("--embedding_type", type=str, default="matrix",
                        choices=["matrix", "vector", "contact"],
                        help="llm embedding type.")
    parser.add_argument("--trunc_type", type=str, default="right",
                        choices=["left", "right"],
                        help="llm trunc type of seq.")
    parser.add_argument("--truncation_seq_length", type=int,
                        default=4094,
                        help="truncation seq length.")
    parser.add_argument("--matrix_add_special_token", action="store_true",
                        help="whether to add special token embedding in seq representation matrix")

    parser.add_argument("--embedding_complete",  action="store_true",
                        help="when the seq len > inference_max_len, then the embedding matrix is completed by segment")
    parser.add_argument("--embedding_complete_seg_overlap",  action="store_true",
                        help="segment overlap")
    parser.add_argument("--embedding_fixed_len_a_time", type=int, default=None,
                        help="the embedding fixed length of once inference for longer sequence")

    parser.add_argument('--gpu', type=int, default=-1, help="the gpu id to use.")

    input_args = parser.parse_args()
    return input_args


if __name__ == "__main__":
    run_args = get_args()
    main(run_args)