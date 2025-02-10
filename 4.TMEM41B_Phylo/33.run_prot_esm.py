# nohup python run_prot_esm.py --seq_type prot > run_prot_esm_all.out 2>&1 &


from inference_embedding_esm import main, get_args

def setup_args():
    args = get_args()
    args.llm_type = "esm2"
    args.llm_version = "3B"

    args.input_file = "/home/share/huadjyin/home/zhoumin3/LucaOne/LucaOneApp-master/data/new/prot_nox.csv"
    args.id_idx = 4
    args.seq_idx = 3
    args.save_path = "/home/share/huadjyin/home/zhoumin3/LucaOne/LucaOneApp-master/Output/new/esm/"

    args.seq_type = "prot"
    args.embedding_type = "matrix"
    args.matrix_add_special_token = True
    args.embedding_complete = True
    args.embedding_complete_seg_overlap = True

    args.truncation_seq_length = 100000
    args.trunc_type = "right"
    args.gpu_id = 1
    args.vector_type = "mean"

    return args

if __name__ == "__main__":
    args = setup_args()

    args_dict = vars(args)
    print("-" * 20 + "Input Args:" + "-" * 20)
    for key, value in sorted(args_dict.items()):
        print(f"{key}: {value}")
    print("-" * 50)

    main(args)
