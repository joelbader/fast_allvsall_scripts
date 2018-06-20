#!/usr/bin/python3

import pickle
import pandas as pd
import xlsxwriter

#Parameters and Settings:
req_match_fraction = 0.75 #Any matches below this fraction identity will not be copied from the pickle to the output excel file
pickle_filename = "mhdb_50.pickle" #Input *.pickle file
output_excel_filename = 'mhdb_'+str(round(100*req_match_fraction))+'.xlsx'

def main():
    with open(pickle_filename, 'rb') as fp:
        db = pickle.load(fp)

    #Write out formatted excel table
    workbook = xlsxwriter.Workbook(output_excel_filename)
    worksheet = workbook.add_worksheet()

    color_ls=['blue','cyan','gray','lime','magenta','orange','pink','red','yellow']
    color_ls=['#3333FF','cyan','silver','lime','#CC00CC','#FF9900','#FF99FF','#FF3333','yellow']
    #Defining Format objects into a list is a workaround to allow for multiple formats until Xlsxwriter developers fix this
    first_col_format_ls = list()
    for i,c in enumerate(db.columns):
        new_format = workbook.add_format()
        new_format.set_bold()
        new_format.set_bg_color(color_ls[i % len(color_ls)])
        first_col_format_ls.append(new_format)

    first_row_format_ls = list()
    for i,r in enumerate(db.index):
        new_format = workbook.add_format()
        new_format.set_bold()
        new_format.set_bg_color(color_ls[i % len(color_ls)])
        first_row_format_ls.append(new_format)

    #Print homologs
    xrow = 0
    xcol = 2
    restart_row = 0 #This is used as the shift for the table rows when moving query genomes
    cols = db.columns
    for i,c in enumerate(cols):
        worksheet.write(xrow, xcol, c[8:-3], first_row_format_ls[i])
        xcol += 1
    xrow += 1
    restart_row=xrow
    for i,r in enumerate(db.index):
        xcol = 0
        #Print first two columns of final table
        proto = db.loc[r,cols[0]]
        for key in proto:
            worksheet.write(xrow, xcol, r[8:-3], first_col_format_ls[i])
            worksheet.write(xrow, xcol+1, key, first_col_format_ls[i])
            xrow += 1
        xcol += 2
        #Print the other columns of the final table
        for c in cols:
            d = db.loc[r,c]
            xrow = restart_row
            for key in proto:
                matches = []
                for match in d[key]:
                    if match[1] > req_match_fraction:
                        matches.append(match[0]+' ('+str(round(match[1],2))+')')
                worksheet.write(xrow, xcol, '\n'.join(matches))
                xrow += 1
            xcol += 1
        restart_row=xrow

    workbook.close()

if __name__=='__main__':
    main()
