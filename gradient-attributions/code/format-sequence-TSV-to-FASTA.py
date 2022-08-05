# %%
import argparse

# %%
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('tsv')
    args = parser.parse_args()

    with open(args.tsv) as f:
        for line in f:
            row = line.strip().split()
            print(f'>{row[0]}:{row[2]}')
            print(row[1])

# %%
if __name__ == '__main__':
    main()