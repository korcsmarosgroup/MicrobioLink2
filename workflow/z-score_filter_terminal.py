import argparse
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde

# def log2_converter(num):
#     try:
#         num = float(num)
#     except:
#         return 'NaN'
#     if num == 0:
#         return 'NaN'
#     elif num > 0:
#         result = np.log2(num)
#         return result

def main(args):
    # Open nested dictionary as a pandas df
    df = pd.read_csv(args.input_file, index_col=0)

    # Save columns (cell-condition) and rows (genes)
    columns = df.columns
    rows = df.index.values

    # Empty df for the log2 TPM values
    #df_cells_log2 = pd.DataFrame(columns=columns, index=rows)

    # Filling the table
    # for column in columns:
    #     df_cells_log2[column] = df[column].apply(log2_converter)

    # Empty df for filtered results
    df_cells_log2_filtered = pd.DataFrame(columns=columns, index=rows)

    # Going through each cell-condition column and apply z-score filter
    for i in columns:
        count = df[i].tolist()
        count = np.array(count, dtype=float)  # Explicitly convert to NumPy array with dtype=float

        count_filtered = count[np.logical_not(np.isnan(count))]

        # creating the Gauss-curve
        kernel = gaussian_kde(count_filtered)

        # creating the X axis -> divide the list for 100 units, xi numpy lists contains that 100 values
        # expected value = most often value, middle of Gaus-curve
        xi = np.linspace(count_filtered.min(), count_filtered.max(), 100)

        # calculate y for each x point
        yi = kernel.evaluate(xi)

        # expected value calculation, which x is by the max y value? (np.argmax(yi) = position)
        mu = xi[np.argmax(yi)]

        # count > mu  = list of boolean values; mean of values right from the expected value
        U = count_filtered[count_filtered > mu].mean()

        # calculation of standard deviation
        sigma = (U - mu) * np.sqrt(np.pi / 2)

        # new score: deviancy from the mean divided by sigma (standard deviation)
        # z-value: relative value - deviation from the mean in the st.dev in the data
        zcount = (count - mu) / sigma

        score_list = [count[list(zcount).index(x)] if x > args.zscore else 'NaN' for x in zcount]

        s = pd.Series(score_list, index=rows)

        df_cells_log2_filtered[i] = s

    # Writing out results
    df_cells_log2_filtered.to_csv(args.output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gene expression filtration based on individual cell count and z-score.")
    parser.add_argument("-i", "--input_file", required=True, help="Input CSV file with gene expression data.")
    parser.add_argument("-zscore", "--zscore", required=True, default = -3, help="Z-score cut-off to filter lowely expressed genes")
    parser.add_argument("-o", "--output_file", required=True, help="Output CSV file for filtered results.")
    args = parser.parse_args()
    main(args)
