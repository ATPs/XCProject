# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 17:20:22 2016

@author: k
"""
def mongoExcel():
    folder = "E:\\Study\\coursera\\Udacity\\ud032\\"
    filename = folder + '2013_ERCOT_Hourly_Load_Data.xls'
    
    import xlrd
    from zipfile import ZipFile
    datafile = filename
    
    
    def open_zip(datafile):
        with ZipFile('{0}.zip'.format(datafile), 'r') as myzip:
            myzip.extractall()
    
    
    def parse_file(datafile):
        workbook = xlrd.open_workbook(datafile)
        sheet = workbook.sheet_by_index(0)
    
        ### example on how you can get the data
        #sheet_data = [[sheet.cell_value(r, col) for col in range(sheet.ncols)] for r in range(sheet.nrows)]
    
        ### other useful methods:
        # print "\nROWS, COLUMNS, and CELLS:"
        # print "Number of rows in the sheet:", 
        # print sheet.nrows
        # print "Type of data in cell (row 3, col 2):", 
        # print sheet.cell_type(3, 2)
        # print "Value in cell (row 3, col 2):", 
        # print sheet.cell_value(3, 2)
        # print "Get a slice of values in column 3, from rows 1-3:"
        # print sheet.col_values(3, start_rowx=1, end_rowx=4)
    
        # print "\nDATES:"
        # print "Type of data in cell (row 1, col 0):", 
        # print sheet.cell_type(1, 0)
        # exceltime = sheet.cell_value(1, 0)
        # print "Time in Excel format:",
        # print exceltime
        # print "Convert time to a Python datetime tuple, from the Excel float:",
        # print xlrd.xldate_as_tuple(exceltime, 0)
        sheet_data = [[sheet.cell_value(r, col) for col in range(sheet.ncols)] for r in range(sheet.nrows)]
        coastdata = [float(i[1]) for i in sheet_data[1:]]
        import numpy as np
        coastdata = np.array(coastdata)
        timedata = [xlrd.xldate_as_tuple(i[0],0) for i in sheet_data[1:]]
        idmax = coastdata.argmax()
        idmin = coastdata.argmin()
        meancoast = coastdata.mean()
        data = {
                'maxtime': timedata[idmax],
                'maxvalue': coastdata[idmax],
                'mintime': timedata[idmin],
                'minvalue': coastdata[idmin],
                'avgcoast': meancoast
        }
        return data
    
    
    def test():
        open_zip(datafile)
        data = parse_file(datafile)
    
        assert data['maxtime'] == (2013, 8, 13, 17, 0, 0)
        assert round(data['maxvalue'], 10) == round(18779.02551, 10)
    
    
    test()

def mongoJson():
    import json
    import requests
    
    
    BASE_URL = "http://musicbrainz.org/ws/2/"
    ARTIST_URL = BASE_URL + "artist/"
    
    # query parameters are given to the requests.get function as a dictionary; this
    # variable contains some starter parameters.
    query_type = {  "simple": {},
                    "atr": {"inc": "aliases+tags+ratings"},
                    "aliases": {"inc": "aliases"},
                    "releases": {"inc": "releases"}}
    
    
    def query_site(url, params, uid="", fmt="json"):
        # This is the main function for making queries to the musicbrainz API.
        # A json document should be returned by the query.
        params["fmt"] = fmt
        r = requests.get(url + uid, params=params)
        print "requesting", r.url
    
        if r.status_code == requests.codes.ok:
            return r.json()
        else:
            r.raise_for_status()
    
    
    def query_by_name(url, params, name):
        # This adds an artist name to the query parameters before making
        # an API call to the function above.
        params["query"] = "artist:" + name
        return query_site(url, params)
    
    
    def pretty_print(data, indent=4):
        # After we get our output, we can format it to be more readable
        # by using this function.
        if type(data) == dict:
            print json.dumps(data, indent=indent, sort_keys=True)
        else:
            print data
    
    
    def main():
        '''
        Modify the function calls and indexing below to answer the questions on
        the next quiz. HINT: Note how the output we get from the site is a
        multi-level JSON document, so try making print statements to step through
        the structure one level at a time or copy the output to a separate output
        file.
        '''
        results = query_by_name(ARTIST_URL, query_type["simple"], "Nirvana")
        pretty_print(results)
    
        artist_id = results["artists"][1]["id"]
        print "\nARTIST:"
        pretty_print(results["artists"][1])
    
        artist_data = query_site(ARTIST_URL, query_type["releases"], artist_id)
        releases = artist_data["releases"]
        print "\nONE RELEASE:"
        pretty_print(releases[0], indent=2)
        release_titles = [r["title"] for r in releases]
    
        print "\nALL TITLES:"
        for t in release_titles:
            print t
    
    
    if __name__ == '__main__':
        main()


def useCSV():
    #!/usr/bin/env python
    """
    Your task is to process the supplied file and use the csv module to extract data from it.
    The data comes from NREL (National Renewable Energy Laboratory) website. Each file
    contains information from one meteorological station, in particular - about amount of
    solar and wind energy for each hour of day.
    
    Note that the first line of the datafile is neither data entry, nor header. It is a line
    describing the data source. You should extract the name of the station from it.
    
    The data should be returned as a list of lists (not dictionaries).
    You can use the csv modules "reader" method to get data in such format.
    Another useful method is next() - to get the next line from the iterator.
    You should only change the parse_file function.
    """
    import csv
    import os
    
    DATADIR = ""
    DATAFILE = "745090.csv"
    
    
    def parse_file(datafile):
        name = ""
        data = []
        with open(datafile,'rb') as f:
    
            a = csv.reader(f)
            name = a.next()[1]
            a.next()
            for row in a:
                data.append(row)
        # Do not change the line below
        return (name, data)
    
    
    def test():
        datafile = os.path.join(DATADIR, DATAFILE)
        name, data = parse_file(datafile)
    
        assert name == "MOUNTAIN VIEW MOFFETT FLD NAS"
        assert data[0][1] == "01:00"
        assert data[2][0] == "01/01/2005"
        assert data[2][5] == "2"
    
    
    if __name__ == "__main__":
        test()

def excel_csv():
    # -*- coding: utf-8 -*-
    '''
    Find the time and value of max load for each of the regions
    COAST, EAST, FAR_WEST, NORTH, NORTH_C, SOUTHERN, SOUTH_C, WEST
    and write the result out in a csv file, using pipe character | as the delimiter.
    
    An example output can be seen in the "example.csv" file.
    '''
    
    import xlrd
    import os
    import csv
    from zipfile import ZipFile
    
    datafile = "2013_ERCOT_Hourly_Load_Data.xls"
    outfile = "2013_Max_Loads.csv"
    
    
    def open_zip(datafile):
        with ZipFile('{0}.zip'.format(datafile), 'r') as myzip:
            myzip.extractall()
    
    
    def parse_file(datafile):
        workbook = xlrd.open_workbook(datafile)
        sheet = workbook.sheet_by_index(0)
        sheet_data = [[sheet.cell_value(r,col) for col in range(sheet.ncols)] for r in range(sheet.nrows)]
        import pandas as pd
        data =[]
        dfsheet = pd.DataFrame(sheet_data[1:],columns=sheet_data[0])
        for column in dfsheet.columns[1:-1]:
            num = dfsheet[[column]].idxmax()
            times = dfsheet[[dfsheet.columns[0]]].values
            time = float(times[num])
            time = xlrd.xldate_as_tuple(time,0)
            data.append({})
            data[-1][column] = {}
            data[-1][column]['Max Load'] = float(dfsheet[[column]].values[num])
            data[-1][column]['Year'] = time[0]
            data[-1][column]['Month'] = time[1]
            data[-1][column]['Day'] = time[2]
            data[-1][column]['Hour'] = time[3]
            data[-1][column]['Minute'] = time[4]
            data[-1][column]['Second'] = time[5]
            data[-1][column]['Station'] = column
            #data.append([column]+ list(time) + [float(dfsheet[[column]].values[num])])
        # YOUR CODE HERE
        # Remember that you can use xlrd.xldate_as_tuple(sometime, 0) to convert
        # Excel date to Python tuple of (year, month, day, hour, minute, second)
        return data
    
    def save_file(data, filename):
        # YOUR CODE HERE
        
        with open(filename,'wb') as f:
            w = csv.DictWriter(f,fieldnames=['Station','Year','Month','Day','Hour','Minute','Second','Max Load'],delimiter='|')
            w.writeheader()
            for ele in data:
                w.writerow(ele.values()[0])
                
            
    
        
    def test():
        open_zip(datafile)
        data = parse_file(datafile)
        save_file(data, outfile)
    
        number_of_rows = 0
        stations = []
    
        ans = {'FAR_WEST': {'Max Load': '2281.2722140000024',
                            'Year': '2013',
                            'Month': '6',
                            'Day': '26',
                            'Hour': '17'}}
        correct_stations = ['COAST', 'EAST', 'FAR_WEST', 'NORTH',
                            'NORTH_C', 'SOUTHERN', 'SOUTH_C', 'WEST']
        fields = ['Year', 'Month', 'Day', 'Hour', 'Max Load']
    
        with open(outfile) as of:
            csvfile = csv.DictReader(of, delimiter="|")
            for line in csvfile:
                station = line['Station']
                if station == 'FAR_WEST':
                    for field in fields:
                        # Check if 'Max Load' is within .1 of answer
                        if field == 'Max Load':
                            max_answer = round(float(ans[station][field]), 1)
                            max_line = round(float(line[field]), 1)
                            assert max_answer == max_line
    
                        # Otherwise check for equality
                        else:
                            assert ans[station][field] == line[field]
    
                number_of_rows += 1
                stations.append(station)
    
            # Output should be 8 lines not including header
            assert number_of_rows == 8
    
            # Check Station Names
            assert set(stations) == set(correct_stations)
    
            
    if __name__ == "__main__":
        test()




def wranglingJSON():
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    """
    This exercise shows some important concepts that you should be aware about:
    - using codecs module to write unicode files
    - using authentication with web APIs
    - using offset when accessing web APIs
    
    To run this code locally you have to register at the NYTimes developer site 
    and get your own API key. You will be able to complete this exercise in our UI
    without doing so, as we have provided a sample result.
    
    Your task is to process the saved file that represents the most popular
    articles (by view count) from the last day, and return the following data:
    - list of dictionaries, where the dictionary key is "section" and value is "title"
    - list of URLs for all media entries with "format": "Standard Thumbnail"
    
    All your changes should be in the article_overview function.
    The rest of functions are provided for your convenience, if you want to access
    the API by yourself.
    """
    import json
    import codecs
    import requests
    
    URL_MAIN = "http://api.nytimes.com/svc/"
    URL_POPULAR = URL_MAIN + "mostpopular/v2/"
    API_KEY = { "popular": "6b148d68e2f572c5170ef9f3e696fe9c:15:75147514",
                "article": "fb46817f7cd16f6abfe8c71b24b33eb7:5:75147514"}
    
    
    def get_from_file(kind, period):
        filename = "popular-{0}-{1}.json".format(kind, period)
        with open(filename, "r") as f:
            return json.loads(f.read())
    
    
    def article_overview(kind, period):
        data = get_from_file(kind, period)
        titles = []
        urls =[]
        # YOUR CODE HERE
    
        return (titles, urls)
    
    
    def query_site(url, target, offset):
        # This will set up the query with the API key and offset
        # Web services often use offset paramter to return data in small chunks
        # NYTimes returns 20 articles per request, if you want the next 20
        # You have to provide the offset parameter
        if API_KEY["popular"] == "" or API_KEY["article"] == "":
            print "You need to register for NYTimes Developer account to run this program."
            print "See Intructor notes for information"
            return False
        params = {"api-key": API_KEY[target], "offset": offset}
        r = requests.get(url, params = params)
    
        if r.status_code == requests.codes.ok:
            return r.json()
        else:
            r.raise_for_status()
    
    
    def get_popular(url, kind, days, section="all-sections", offset=0):
        # This function will construct the query according to the requirements of the site
        # and return the data, or print an error message if called incorrectly
        if days not in [1,7,30]:
            print "Time period can be 1,7, 30 days only"
            return False
        if kind not in ["viewed", "shared", "emailed"]:
            print "kind can be only one of viewed/shared/emailed"
            return False
    
        url += "most{0}/{1}/{2}.json".format(kind, section, days)
        data = query_site(url, "popular", offset)
    
        return data
    
    
    def save_file(kind, period):
        # This will process all results, by calling the API repeatedly with supplied offset value,
        # combine the data and then write all results in a file.
        data = get_popular(URL_POPULAR, "viewed", 1)
        num_results = data["num_results"]
        full_data = []
        with codecs.open("popular-{0}-{1}.json".format(kind, period), encoding='utf-8', mode='w') as v:
            for offset in range(0, num_results, 20):        
                data = get_popular(URL_POPULAR, kind, period, offset=offset)
                full_data += data["results"]
            
            v.write(json.dumps(full_data, indent=2))
    
    
    def test():
        titles, urls = article_overview("viewed", 1)
        assert len(titles) == 20
        assert len(urls) == 30
        assert titles[2] == {'Opinion': 'Professors, We Need You!'}
        assert urls[20] == 'http://graphics8.nytimes.com/images/2014/02/17/sports/ICEDANCE/ICEDANCE-thumbStandard.jpg'
    
    
    if __name__ == "__main__":
        test()

def xmlExtractingData():
    #!/usr/bin/env python
    # Your task here is to extract data from xml on authors of an article
    # and add it to a list, one item for an author.
    # See the provided data structure for the expected format.
    # The tags for first name, surname and email should map directly
    # to the dictionary keys
    import xml.etree.ElementTree as ET
    
    article_file = r"E:\Study\coursera\Udacity\ud032\Lesson_2_Data_in_More_Complex_Formats\07-Extracting_Data\exampleResearchArticle.xml"
    
    
    def get_root(fname):
        tree = ET.parse(fname)
        return tree.getroot()
    
    
    def get_authors(root):
        authors = []
        for author in root.findall('./fm/bibl/aug/au'):
            data = {
                    "fnm": None,
                    "snm": None,
                    "email": None
            }
            data['fnm'] = author.find('./fnm').text
            data['snm'] = author.find('./snm').text
            data['email'] = author.find('./email').text
            iids = author.findall('./insr')
            for iid in iids:
                data['insr'].append(iid.attrib.values()[0])
            # YOUR CODE HERE
    
            authors.append(data)
    
        return authors
    
    
    def test():
        solution = [{'fnm': 'Omer', 'snm': 'Mei-Dan', 'email': 'omer@extremegate.com'}, {'fnm': 'Mike', 'snm': 'Carmont', 'email': 'mcarmont@hotmail.com'}, {'fnm': 'Lior', 'snm': 'Laver', 'email': 'laver17@gmail.com'}, {'fnm': 'Meir', 'snm': 'Nyska', 'email': 'nyska@internet-zahav.net'}, {'fnm': 'Hagay', 'snm': 'Kammar', 'email': 'kammarh@gmail.com'}, {'fnm': 'Gideon', 'snm': 'Mann', 'email': 'gideon.mann.md@gmail.com'}, {'fnm': 'Barnaby', 'snm': 'Clarck', 'email': 'barns.nz@gmail.com'}, {'fnm': 'Eugene', 'snm': 'Kots', 'email': 'eukots@gmail.com'}]
        
        root = get_root(article_file)
        data = get_authors(root)
    
        assert data[0] == solution[0]
        assert data[1]["fnm"] == solution[1]["fnm"]
    
    
    test()




def usingBeautifulSoup():
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    # Please note that the function 'make_request' is provided for your reference only.
    # You will not be able to to actually use it from within the Udacity web UI.
    # Your task is to process the HTML using BeautifulSoup, extract the hidden
    # form field values for "__EVENTVALIDATION" and "__VIEWSTATE" and set the appropriate
    # values in the data dictionary.
    # All your changes should be in the 'extract_data' function
    from bs4 import BeautifulSoup
    import requests
    import json
    
    html_page = r"E:\Study\coursera\Udacity\ud032\Lesson_2_Data_in_More_Complex_Formats\18-Using_Beautiful_Soup\page_source.html"
    
    
    def extract_data(page):
        data = {"eventvalidation": "",
                "viewstate": ""}
        with open(page, "r") as html:
            # do something here to find the necessary values
            soup = BeautifulSoup(html,'html.parser')
            ev = soup.find(id="__EVENTVALIDATION")
            data["eventvalidation"] = ev["value"]

            vs = soup.find(id="__VIEWSTATE")
            data["viewstate"] = vs["value"]
                
    
        return data
    
    
    def make_request(data):
        eventvalidation = data["eventvalidation"]
        viewstate = data["viewstate"]
    
        r = requests.post("http://www.transtats.bts.gov/Data_Elements.aspx?Data=2",
                        data={'AirportList': "BOS",
                              'CarrierList': "VX",
                              'Submit': 'Submit',
                              "__EVENTTARGET": "",
                              "__EVENTARGUMENT": "",
                              "__EVENTVALIDATION": eventvalidation,
                              "__VIEWSTATE": viewstate
                        })
    
        return r.text
    
    
    def test():
        data = extract_data(html_page)
        assert data["eventvalidation"] != ""
        assert data["eventvalidation"].startswith("/wEWjAkCoIj1ng0")
        assert data["viewstate"].startswith("/wEPDwUKLTI")
    
        
    test()


def carrierList():
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    """
    Please note that the function 'make_request' is provided for your reference only.
    You will not be able to to actually use it from within the Udacity web UI.
    All your changes should be in the 'extract_carrier' function.
    Also note that the html file is a stripped down version of what is actually on
    the website.
    
    Your task in this exercise is to get a list of all airlines. Exclude all of the
    combination values like "All U.S. Carriers" from the data that you return.
    You should return a list of codes for the carriers.
    """
    
    from bs4 import BeautifulSoup
    import requests
    html_page = r"E:\Study\coursera\Udacity\ud032\Lesson_2_Problem_Set\01-Carrier_List\options.html"
    
    
    def extract_carriers(page):
        data = []
    
        with open(page, "r") as html:
            # do something here to find the necessary values
            soup = BeautifulSoup(html, "lxml")
            carrierlist = soup.find(id = 'CarrierList')
            carriers = carrierlist.findAll('option')
            for carrier in carriers:
                data.append(carrier['value'])
    
        return data
    
    
    def extract_airports(page):
        data = []
        with open(page, "r") as html:
            # do something here to find the necessary values
            soup = BeautifulSoup(html, "lxml")
            airportlist = soup.find(id='AirportList')
            airports = airportlist.findAll("option")
            for airport in airports:
                data.append(airport['value'])
            
    
        return data[3:-3] + data[-2:]
    
    def make_request(data):
        eventvalidation = data["eventvalidation"]
        viewstate = data["viewstate"]
        airport = data["airport"]
        carrier = data["carrier"]
    
        r = requests.post("http://www.transtats.bts.gov/Data_Elements.aspx?Data=2",
                        data={'AirportList': airport,
                              'CarrierList': carrier,
                              'Submit': 'Submit',
                              "__EVENTTARGET": "",
                              "__EVENTARGUMENT": "",
                              "__EVENTVALIDATION": eventvalidation,
                              "__VIEWSTATE": viewstate
                        })
    
        return r.text
    
    
    def test():
        data = extract_carriers(html_page)
        assert len(data) == 16
        assert "FL" in data
        assert "NK" in data
    
    test()


def processingAll():
    
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    """
    Let's assume that you combined the code from the previous 2 exercises with code
    from the lesson on how to build requests, and downloaded all the data locally.
    The files are in a directory "data", named after the carrier and airport:
    "{}-{}.html".format(carrier, airport), for example "FL-ATL.html".
    
    The table with flight info has a table class="dataTDRight". Your task is to
    extract the flight data from that table as a list of dictionaries, each
    dictionary containing relevant data from the file and table row. This is an
    example of the data structure you should return:
    
    data = [{"courier": "FL",
             "airport": "ATL",
             "year": 2012,
             "month": 12,
             "flights": {"domestic": 100,
                         "international": 100}
            },
             {"courier": "..."}
    ]
    
    Note - year, month, and the flight data should be integers.
    You should skip the rows that contain the TOTAL data for a year.
    
    There are couple of helper functions to deal with the data files.
    Please do not change them for grading purposes.
    All your changes should be in the 'process_file' function.
    """
    from bs4 import BeautifulSoup
    from zipfile import ZipFile
    import os
    
    datadir = "data"
    
    
    def open_zip(datadir):
        with ZipFile('{0}.zip'.format(datadir), 'r') as myzip:
            myzip.extractall()
    
    
    def process_all(datadir):
        files = os.listdir(datadir)
        return files
    
    
    def process_file(f):
        """
        This function extracts data from the file given as the function argument in
        a list of dictionaries. This is example of the data structure you should
        return:
    
        data = [{"courier": "FL",
                 "airport": "ATL",
                 "year": 2012,
                 "month": 12,
                 "flights": {"domestic": 100,
                             "international": 100}
                },
                {"courier": "..."}
        ]
    
    
        Note - year, month, and the flight data should be integers.
        You should skip the rows that contain the TOTAL data for a year.
        """
        data = []
        info = {}
        info["courier"], info["airport"] = f[:6].split("-")
        # Note: create a new dictionary for each entry in the output data list.
        # If you use the info dictionary defined here each element in the list 
        # will be a reference to the same info dictionary.
        with open("{}/{}".format(datadir, f), "r") as html:
    
            soup = BeautifulSoup(html)
            table = soup.find(id='DataGrid1')
            rows = table.find_all('tr',class_='dataTDRight')
            for row in rows:
                eles = row.findAll('td')
                if eles[1].text != "TOTAL":
                    info['year'] = int(eles[0].text)
                    info['month'] = int(eles[1].text)
                    info['flights'] ={}
                    info['flights']['domestic'] = int(eles[2].text.replace(",",""))
                    info['flights']['international'] = int(eles[3].text.replace(",",""))
                    data.append(info)
                    info={}
                    info["courier"], info["airport"] = f[:6].split("-")
    
        return data
    
    
    def test():
        print "Running a simple test..."
        open_zip(datadir)
        files = process_all(datadir)
        data = []
        # Test will loop over three data files.
        for f in files:
            data += process_file(f)
            
        assert len(data) == 399  # Total number of rows
        for entry in data[:3]:
            assert type(entry["year"]) == int
            assert type(entry["month"]) == int
            assert type(entry["flights"]["domestic"]) == int
            assert len(entry["airport"]) == 3
            assert len(entry["courier"]) == 2
        assert data[0]["courier"] == 'FL'
        assert data[0]["month"] == 10
        assert data[-1]["airport"] == "ATL"
        assert data[-1]["flights"] == {'international': 108289, 'domestic': 701425}
        
        print "... success!"
    
    if __name__ == "__main__":
    test()




def patent_database():
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    """
    This and the following exercise are using US Patent database. The patent.data
    file is a small excerpt of much larger datafiles that are available for
    download from US Patent website. These files are pretty large ( >100 MB each).
    The original file is ~600MB large, you might not be able to open it in a text
    editor.
    
    The data itself is in XML, however there is a problem with how it's formatted.
    Please run this script and observe the error. Then find the line that is
    causing the error. You can do that by just looking at the datafile in the web
    UI, or programmatically. For quiz purposes it does not matter, but as an
    exercise we suggest that you try to do it programmatically.
    
    NOTE: You do not need to correct the error - for now, just find where the error
    is occurring.
    """
    
    import xml.etree.ElementTree as ET
    
    PATENTS = r'E:\Study\coursera\Udacity\ud032\Lesson_2_Problem_Set\04-Patent_Database\patent.data'
    
    def get_root(fname):
    
        tree = ET.parse(fname)
        return tree.getroot()
    
    
    get_root(PATENTS)
    
    
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    # So, the problem is that the gigantic file is actually not a valid XML, because
    # it has several root elements, and XML declarations.
    # It is, a matter of fact, a collection of a lot of concatenated XML documents.
    # So, one solution would be to split the file into separate documents,
    # so that you can process the resulting files as valid XML documents.
    
    import xml.etree.ElementTree as ET
    PATENTS = r'E:\Study\coursera\Udacity\ud032\Lesson_2_Problem_Set\04-Patent_Database\patent.data'
    
    def get_root(fname):
        tree = ET.parse(fname)
        return tree.getroot()
    
    
    def split_file(filename):
        """
        Split the input file into separate files, each containing a single patent.
        As a hint - each patent declaration starts with the same line that was
        causing the error found in the previous exercises.
        
        The new files should be saved with filename in the following format:
        "{}-{}".format(filename, n) where n is a counter, starting from 0.
        """
        fo = open(filename,'r')
        content = fo.read()
        fo.close()
        ls_content = content.split("<?xml")[1:]
        for num in range(len(ls_content)):
            fout = open(filename+str(num),'w')
            fout.write("<?xml"+ls_content[num])
            fout.close()
        

    
    
    def test():
        split_file(PATENTS)
        for n in range(4):
            try:
                fname = "{}-{}".format(PATENTS, n)
                f = open(fname, "r")
                if not f.readline().startswith("<?xml"):
                    print "You have not split the file {} in the correct boundary!".format(fname)
                f.close()
            except:
                print "Could not find file {}. Check if the filename is correct!".format(fname)
    
    
    test()

def correctingValidity():
    """
    Your task is to check the "productionStartYear" of the DBPedia autos datafile for valid values.
    The following things should be done:
    - check if the field "productionStartYear" contains a year
    - check if the year is in range 1886-2014
    - convert the value of the field to be just a year (not full datetime)
    - the rest of the fields and values should stay the same
    - if the value of the field is a valid year in the range as described above,
      write that line to the output_good file
    - if the value of the field is not a valid year as described above, 
      write that line to the output_bad file
    - discard rows (neither write to good nor bad) if the URI is not from dbpedia.org
    - you should use the provided way of reading and writing data (DictReader and DictWriter)
      They will take care of dealing with the header.
    
    You can write helper functions for checking the data and writing the files, but we will call only the 
    'process_file' with 3 arguments (inputfile, output_good, output_bad).
    """
    import csv
    import pprint
    
    INPUT_FILE = r'E:\Study\coursera\Udacity\ud032\Lesson_3_Data_Quality\12-Correcting_Validity\autos.csv'
    OUTPUT_GOOD = 'autos-valid.csv'
    OUTPUT_BAD = 'FIXME-autos.csv'
    
    def process_file(input_file, output_good, output_bad):
        
        ls_good=[]
        ls_bad=[]
        import re
    
        with open(input_file, "r") as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames
            
            for row in reader:
                if 'http://dbpedia.org/' in row['URI']:
                    if re.match('^\d\d\d\d-', row['productionStartYear']):
                        year = int(row['productionStartYear'][:4])
                        if year >=1886 and year <= 2014:
                            row['productionStartYear'] = row['productionStartYear'][:4]
                            ls_good.append(row)
                        else:
                            ls_bad.append(row)
                    else:
                        ls_bad.append(row)
    
            #COMPLETE THIS FUNCTION
    
    
    
        # This is just an example on how you can use csv.DictWriter
        # Remember that you have to output 2 files
        with open(output_good, "w") as g:
            writer = csv.DictWriter(g, delimiter=",", fieldnames= header)
            writer.writeheader()
            for row in ls_good:
                writer.writerow(row)
        
        with open(output_bad, "w") as h:
            writer = csv.DictWriter(h, delimiter=",", fieldnames= header)
            writer.writeheader()
            for row in ls_bad:
                writer.writerow(row)
    
    
    def test():
    
        process_file(INPUT_FILE, OUTPUT_GOOD, OUTPUT_BAD)
    
    
    if __name__ == "__main__":
        test()


def auditingDataQuality():
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    """
    In this problem set you work with cities infobox data, audit it, come up with a
    cleaning idea and then clean it up. In the first exercise we want you to audit
    the datatypes that can be found in some particular fields in the dataset.
    The possible types of values can be:
    - NoneType if the value is a string "NULL" or an empty string ""
    - list, if the value starts with "{"
    - int, if the value can be cast to int
    - float, if the value can be cast to float, but CANNOT be cast to int.
       For example, '3.23e+07' should be considered a float because it can be cast
       as float but int('3.23e+07') will throw a ValueError
    - 'str', for all other values
    
    The audit_file function should return a dictionary containing fieldnames and a 
    SET of the types that can be found in the field. e.g.
    {"field1": set([type(float()), type(int()), type(str())]),
     "field2": set([type(str())]),
      ....
    }
    The type() function returns a type object describing the argument given to the 
    function. You can also use examples of objects to create type objects, e.g.
    type(1.1) for a float: see the test function below for examples.
    
    Note that the first three rows (after the header row) in the cities.csv file
    are not actual data points. The contents of these rows should note be included
    when processing data types. Be sure to include functionality in your code to
    skip over or detect these rows.
    """
    import codecs
    import csv
    import json
    import pprint
    
    CITIES = r'E:\Study\coursera\Udacity\ud032\Lesson_3_Problem_Set\01-Auditing_Data_Quality\cities.csv'
    
    FIELDS = ["name", "timeZone_label", "utcOffset", "homepage", "governmentType_label",
              "isPartOf_label", "areaCode", "populationTotal", "elevation",
              "maximumElevation", "minimumElevation", "populationDensity",
              "wgs84_pos#lat", "wgs84_pos#long", "areaLand", "areaMetro", "areaUrban"]
    
    def audit_file(filename, fields):
        fieldtypes = {}
        
        for field in fields:
            fieldtypes[field] = set()
        f = open(filename,'r')
        dc_f = csv.DictReader(f)
        ls_f =[row for row in dc_f]
        for ele in ls_f[3:]:
            
            for field in fields:
                e = ele[field]
                if e == 'NULL' or e == '':
                    fieldtypes[field].add(type(None))
                elif e[0] == '{':
                    fieldtypes[field].add(type([]))
                else:
                    try:
                        int(e)
                        fieldtypes[field].add(type(1))
                    except:
                        try:
                            float(e)
                            fieldtypes[field].add(type(1.1))
                        except:
                            fieldtypes[field].add(type('aa'))
    
        # YOUR CODE HERE
    
    
        return fieldtypes
    
    
    def test():
        fieldtypes = audit_file(CITIES, FIELDS)
    
        pprint.pprint(fieldtypes)
    
        assert fieldtypes["areaLand"] == set([type(1.1), type([]), type(None)])
        assert fieldtypes['areaMetro'] == set([type(1.1), type(None)])
        
    if __name__ == "__main__":
        test()



def mongoDBproblemset3():
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    """
    In this problem set you work with another type of infobox data, audit it,
    clean it, come up with a data model, insert it into MongoDB and then run some
    queries against your database. The set contains data about Arachnid class
    animals.
    
    Your task in this exercise is to parse the file, process only the fields that
    are listed in the FIELDS dictionary as keys, and return a list of dictionaries
    of cleaned values. 
    
    The following things should be done:
    - keys of the dictionary changed according to the mapping in FIELDS dictionary
    - trim out redundant description in parenthesis from the 'rdf-schema#label'
      field, like "(spider)"
    - if 'name' is "NULL" or contains non-alphanumeric characters, set it to the
      same value as 'label'.
    - if a value of a field is "NULL", convert it to None
    - if there is a value in 'synonym', it should be converted to an array (list)
      by stripping the "{}" characters and splitting the string on "|". Rest of the
      cleanup is up to you, e.g. removing "*" prefixes etc. If there is a singular
      synonym, the value should still be formatted in a list.
    - strip leading and ending whitespace from all fields, if there is any
    - the output structure should be as follows:
    
    [ { 'label': 'Argiope',
        'uri': 'http://dbpedia.org/resource/Argiope_(spider)',
        'description': 'The genus Argiope includes rather large and spectacular spiders that often ...',
        'name': 'Argiope',
        'synonym': ["One", "Two"],
        'classification': {
                          'family': 'Orb-weaver spider',
                          'class': 'Arachnid',
                          'phylum': 'Arthropod',
                          'order': 'Spider',
                          'kingdom': 'Animal',
                          'genus': None
                          }
      },
      { 'label': ... , }, ...
    ]
    
      * Note that the value associated with the classification key is a dictionary
        with taxonomic labels.
    """
    import codecs
    import csv
    import json
    import pprint
    import re
    
    DATAFILE = r'E:\Study\coursera\Udacity\ud032\Lesson_4_Problem_Set\01-Preparing_Data\arachnid.csv'
    FIELDS ={'rdf-schema#label': 'label',
             'URI': 'uri',
             'rdf-schema#comment': 'description',
             'synonym': 'synonym',
             'name': 'name',
             'family_label': 'family',
             'class_label': 'class',
             'phylum_label': 'phylum',
             'order_label': 'order',
             'kingdom_label': 'kingdom',
             'genus_label': 'genus'}
    
    
    def process_file(filename, fields):
        
        import re
    
        process_fields = fields.keys()
        data = []
        with open(filename, "r") as f:
            reader = csv.DictReader(f)
            for i in range(3):
                l = reader.next()
    
            for line in reader:
                # YOUR CODE HERE
                e = {}
                e['classification'] = {}
                for key, value in line.items():
                    if key in fields:
                        keyn = fields[key]
                    e['uri'] = line["URI"]
                    e['description'] = line['rdf-schema#comment']
                    if key == 'rdf-schema#label':
                        value = value.split(' (')[0]
                        e[keyn] = value.strip()
                    elif key == 'name':
                        if value == 'NULL' or not re.match('[\w ]', value):
                            value = line['rdf-schema#label'].split(' (')[0]
                        e[keyn] = value.strip()
                    else:
                        if value == 'NULL':
                            value = None
                        else:
                            value = value.strip()
                        if key == 'synonym':
                            if value != None and value != "":
                                if value[0] == '{':
                                    e['synonym'] = [i.strip() for i in value[1:-1].split('|')]
                                else:
                                    e['synonym'] = [value.strip()]
                            else:
                                e['synonym'] = None
                        if key == 'family_label':
                            e['classification']['family'] = value
                        if key == 'class_label':
                            e['classification']['class'] = value
                        if key == 'phylum_label':
                            e['classification']['phylum'] = value
                        if key == 'order_label':
                            e['classification']['order'] = value
                        if key == 'kingdom_label':
                            e['classification']['kingdom'] = value
                        if key == 'genus_label':
                            e['classification']['genus'] = value
                    
                
                data.append(e)
                
        return data
    
    
    def parse_array(v):
        if (v[0] == "{") and (v[-1] == "}"):
            v = v.lstrip("{")
            v = v.rstrip("}")
            v_array = v.split("|")
            v_array = [i.strip() for i in v_array]
            return v_array
        return [v]
    
    
    def test():
        data = process_file(DATAFILE, FIELDS)
        print "Your first entry:"
        pprint.pprint(data[0])
        first_entry = {
            "synonym": None, 
            "name": "Argiope", 
            "classification": {
                "kingdom": "Animal", 
                "family": "Orb-weaver spider", 
                "order": "Spider", 
                "phylum": "Arthropod", 
                "genus": None, 
                "class": "Arachnid"
            }, 
            "uri": "http://dbpedia.org/resource/Argiope_(spider)", 
            "label": "Argiope", 
            "description": "The genus Argiope includes rather large and spectacular spiders that often have a strikingly coloured abdomen. These spiders are distributed throughout the world. Most countries in tropical or temperate climates host one or more species that are similar in appearance. The etymology of the name is from a Greek name meaning silver-faced."
        }
        print data[17]["name"],data[48]["label"],data[14]["synonym"]
    
        assert len(data) == 76
        assert data[0] == first_entry
        assert data[17]["name"] == "Ogdenia"
        assert data[48]["label"] == "Hydrachnidiae"
        assert data[14]["synonym"] == ["Cyrene Peckham & Peckham"]
    
    if __name__ == "__main__":
        test()


def analyzyingData():
    #!/usr/bin/env python
    """
    The tweets in our twitter collection have a field called "source". This field describes the application
    that was used to create the tweet. Following the examples for using the $group operator, your task is 
    to modify the 'make-pipeline' function to identify most used applications for creating tweets. 
    As a check on your query, 'web' is listed as the most frequently used application.
    'Ubertwitter' is the second most used. The number of counts should be stored in a field named 'count'
    (see the assertion at the end of the script).
    
    Please modify only the 'make_pipeline' function so that it creates and returns an aggregation pipeline
    that can be passed to the MongoDB aggregate function. As in our examples in this lesson, the aggregation 
    pipeline should be a list of one or more dictionary objects. 
    Please review the lesson examples if you are unsure of the syntax.
    
    Your code will be run against a MongoDB instance that we have provided. 
    If you want to run this code locally on your machine, you have to install MongoDB, 
    download and insert the dataset.
    For instructions related to MongoDB setup and datasets please see Course Materials.
    
    Please note that the dataset you are using here is a smaller version of the twitter dataset 
    used in examples in this lesson. 
    If you attempt some of the same queries that we looked at in the lesson examples,
    your results will be different.
    """
    
    
    def get_db(db_name):
        from pymongo import MongoClient
        client = MongoClient('localhost:27017')
        db = client[db_name]
        return db
    
    def make_pipeline():
        # complete the aggregation pipeline
        pipeline = [{'$group': {"_id":'$source', 'count': {'$sum':1}}}, {'$sort': {'count': -1}}] # the '_id' is always needed.
        #print pipeline
        return pipeline
    
    def tweet_sources(db, pipeline):
        return [doc for doc in db.tweets.aggregate(pipeline)]
    
    if __name__ == '__main__':
        db = get_db('twitter')
        pipeline = make_pipeline()
        result = tweet_sources(db, pipeline)
        import pprint
        pprint.pprint(result[0])
        assert result[0] == {u'count': 868, u'_id': u'web'}


"""
more operations for aggregate
$project, used to reshape the data, pull out only the interested field.
$match, filter the document
$group, operators including $sum, first, last, max, min, avg, push, addToSet(like set, keep unique)
$sort
$skip, skip several input documents.
$limit, only apply to several documents.
$unwind, for value in array fields, broken the document to documents, based on arrary fields.
"""
def highest_ration():
    result = db.tweets.aggregate([{'$match':{'user.friends_count':{'$gt':0},
                                             'user.followers_count': {'$gt':0}}},
                                    {'$project':{'ratio':{'$divide':['$user.followers_count','$user.friends_count']},
                                                 'screen_name':'$user.screen_name'}},
                                    {'$sort':{'ratio':-1}},
                                    {'$limit':1}])
    result2 = [{'$match':{'user.time_zone':'Brasilia', 'user.statuses_count':{'$gt':100}}},\
               {'$project':{'followers':'$user.followers_count', 'screen_name':'$user.screen_name', 'tweets':'$user.statuses_count'}},\
               {'$sort':{'followers': -1}},\
               {'$limit':1}]
    
    result3 = [{'$unwind': '$entities.user_mentions'},\
                {'$group':{'_id':'$user.screen_name', 'count': {'$sum':1}}},\
                {'$sort': {'count': -1}},\
                {'$limit':1}]
    
    result4 = [{'$unwind': '$isPartOf'},\
                {'$match': {'country': 'India'}},
                {'$group': {'_id':'$isPartOf','count':{'$sum':1}}},\
                {'$sort': {'count': -1}},\
                {'$limit':1}]
    
    result5 = [{'$unwind': '$entities.hashtags'},\
                {'$group':{'_id':'$entities.hashtags.text', 'retweet_avg': {'$avg':'$retweet_count'}}},\
                {'$sort':{'retweet_avg':-1}}]
    
    result6 = [{'$unwind':'$entities.hashtags'},\
               {'$group':{'_id': '$user.screen_name', 'unique_hashtags': {"$addToSet": '$entities.hashtags.text'}}},\
               {'$sort':{'_id': -1}} ]
    
    result7 = [\
               {'$group':{'_id': '$user.screen_name', 'tweet_text': {"$push": '$text'}, 'count':{'$sum':1}}},\
               {'$sort':{'count': -1}},\
               {'$limit':5}]
            
    result8_unique_user_mentions = [{'$unwind':'$entities.user_mentions'},\
               {'$group':{'_id': '$user.screen_name', 'mset': {"$addToSet": '$entities.user_mentions.screen_name'}}},\
               {'$unwind':'$mset'},\
               {'$group':{'_id':"$_id","count":{"$sum":1}}},\
               {'$sort':{'count': -1}},\
               {'$limit':10}]
    
    
    result9_IndiaPopulationRegionAverage = [{'$match': {"country" : "India"}},\
                {'$unwind':'$isPartOf'},\
               {'$group':{'_id': '$isPartOf', 'avg_region': {"$avg": '$population'}}},\
               {'$group':{'_id':"Nothing","avg":{"$avg":'$avg_region'}}}
               ]
    
    result10_CountryPopulationRegionAverage = [\
                {'$unwind':'$isPartOf'},\
               {'$group':{'_id': {'region':'$isPartOf','country':'$country'},'avg_region': {"$avg": '$population'}}},\
               {'$group':{'_id':"$_id.country","avgRegionalPopulation":{"$avg":'$avg_region'}}}
               ]
    
    return result
    


#index; increase reading, but slower writing as you have to update index after changing the database
#db,ensureIndex({'tg':1}) , ensure_index
#geospatial indexes, $near operation is needed.

def mostCommonCityName():
    #!/usr/bin/env python
    """
    Use an aggregation query to answer the following question. 
    
    What is the most common city name in our cities collection?
    
    Your first attempt probably identified None as the most frequently occurring
    city name. What that actually means is that there are a number of cities
    without a name field at all. It's strange that such documents would exist in
    this collection and, depending on your situation, might actually warrant
    further cleaning. 
    
    To solve this problem the right way, we should really ignore cities that don't
    have a name specified. As a hint ask yourself what pipeline operator allows us
    to simply filter input? How do we test for the existence of a field?
    
    Please modify only the 'make_pipeline' function so that it creates and returns
    an aggregation pipeline that can be passed to the MongoDB aggregate function.
    As in our examples in this lesson, the aggregation pipeline should be a list of
    one or more dictionary objects. Please review the lesson examples if you are
    unsure of the syntax.
    
    Your code will be run against a MongoDB instance that we have provided. If you
    want to run this code locally on your machine, you have to install MongoDB, 
    download and insert the dataset. For instructions related to MongoDB setup and
    datasets please see Course Materials.
    
    Please note that the dataset you are using here is a different version of the
    cities collection provided in the course materials. If you attempt some of the
    same queries that we look at in the problem set, your results may be different.
    """
    
    def get_db(db_name):
        from pymongo import MongoClient
        client = MongoClient('localhost:27017')
        db = client[db_name]
        return db
    
    def make_pipeline():
        # complete the aggregation pipeline
        pipeline = [ \
        {'$match':{'name':{'$exists':1}}},\
        {'$group':{'_id':'$name','count':{'$sum':1}}}, \
        {'$sort':{'count':-1}},\
        {'$limit':1}
        ]
        return pipeline
    
    def aggregate(db, pipeline):
        return [doc for doc in db.cities.aggregate(pipeline)]
    
    
    if __name__ == '__main__':
        # The following statements will be used to test your code by the grader.
        # Any modifications to the code past this point will not be reflected by
        # the Test Run.
        db = get_db('examples')
        pipeline = make_pipeline()
        result = aggregate(db, pipeline)
        import pprint
        pprint.pprint(result[0])
        assert len(result) == 1
        assert result[0] == {'_id': 'Shahpur', 'count': 6}


def regionCities():
    #!/usr/bin/env python
    """
    Use an aggregation query to answer the following question. 
    
    Which Region in India has the largest number of cities with longitude between
    75 and 80?
    
    Please modify only the 'make_pipeline' function so that it creates and returns
    an aggregation pipeline that can be passed to the MongoDB aggregate function.
    As in our examples in this lesson, the aggregation pipeline should be a list of
    one or more dictionary objects. Please review the lesson examples if you are
    unsure of the syntax.
    
    Your code will be run against a MongoDB instance that we have provided. If you
    want to run this code locally on your machine, you have to install MongoDB,
    download and insert the dataset. For instructions related to MongoDB setup and
    datasets please see Course Materials.
    
    Please note that the dataset you are using here is a different version of the
    cities collection provided in the course materials. If you attempt some of the
    same queries that we look at in the problem set, your results may be different.
    """
    
    def get_db(db_name):
        from pymongo import MongoClient
        client = MongoClient('localhost:27017')
        db = client[db_name]
        return db
    
    def make_pipeline():
        # complete the aggregation pipeline
        pipeline = [ {'$match': {"country" : "India",'lon':{'$gte':75,'$lte':80}}},\
                {'$unwind':'$isPartOf'},\
               {'$group':{'_id': '$isPartOf', 'count': {"$sum": 1}}},\
               {'$sort':{'count':-1}},\
               {'$limit':1}\
               ]
        return pipeline
    
    def aggregate(db, pipeline):
        return [doc for doc in db.cities.aggregate(pipeline)]
    
    if __name__ == '__main__':
        # The following statements will be used to test your code by the grader.
        # Any modifications to the code past this point will not be reflected by
        # the Test Run.
        db = get_db('examples')
        pipeline = make_pipeline()
        result = aggregate(db, pipeline)
        import pprint
        pprint.pprint(result[0])
        assert len(result) == 1
        assert result[0]["_id"] == 'Tamil Nadu'
        assert result[0]["count"] == 424




def iterativeParsing():
    #for envent, elem in ET.iterparse(osm_file)
    import xml.etree.cElementTree as ET
    import pprint
    filename = r'E:\Study\coursera\Udacity\ud032\chicago_illinois.osm'
    fo = open(filename,'r')
    data = ET.iterparse(fo)
    tags = {}
    for tag, elem in data:
        print tag, elem
        break
    def is_street_name(ele):
        return (elem.attrib['k'] == 'addr:street')
    def audit():
        for event, elem in ET.iterparse(osm_file,events = ("start",)):
            if elem.tag == "way":
                for tag in elem.iter("tag"):
                    if is_street_name(tag):
                        audit_street_type(street_types, tag.attrib['v'])
        pprint.pprint(dict(street_types))


def finalproject():
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    import xml.etree.cElementTree as ET
    import pprint
    import re
    import codecs
    import json
    """
    Your task is to wrangle the data and transform the shape of the data
    into the model we mentioned earlier. The output should be a list of dictionaries
    that look like this:
    
    {
    "id": "2406124091",
    "type: "node",
    "visible":"true",
    "created": {
              "version":"2",
              "changeset":"17206049",
              "timestamp":"2013-08-03T16:43:42Z",
              "user":"linuxUser16",
              "uid":"1219059"
            },
    "pos": [41.9757030, -87.6921867],
    "address": {
              "housenumber": "5157",
              "postcode": "60625",
              "street": "North Lincoln Ave"
            },
    "amenity": "restaurant",
    "cuisine": "mexican",
    "name": "La Cabana De Don Luis",
    "phone": "1 (773)-271-5176"
    }
    
    You have to complete the function 'shape_element'.
    We have provided a function that will parse the map file, and call the function with the element
    as an argument. You should return a dictionary, containing the shaped data for that element.
    We have also provided a way to save the data in a file, so that you could use
    mongoimport later on to import the shaped data into MongoDB. 
    
    Note that in this exercise we do not use the 'update street name' procedures
    you worked on in the previous exercise. If you are using this code in your final
    project, you are strongly encouraged to use the code from previous exercise to 
    update the street names before you save them to JSON. 
    
    In particular the following things should be done:
    - you should process only 2 types of top level tags: "node" and "way"
    - all attributes of "node" and "way" should be turned into regular key/value pairs, except:
        - attributes in the CREATED array should be added under a key "created"
        - attributes for latitude and longitude should be added to a "pos" array,
          for use in geospacial indexing. Make sure the values inside "pos" array are floats
          and not strings. 
    - if the second level tag "k" value contains problematic characters, it should be ignored
    - if the second level tag "k" value starts with "addr:", it should be added to a dictionary "address"
    - if the second level tag "k" value does not start with "addr:", but contains ":", you can
      process it in a way that you feel is best. For example, you might split it into a two-level
      dictionary like with "addr:", or otherwise convert the ":" to create a valid key.
    - if there is a second ":" that separates the type/direction of a street,
      the tag should be ignored, for example:
    
    <tag k="addr:housenumber" v="5158"/>
    <tag k="addr:street" v="North Lincoln Avenue"/>
    <tag k="addr:street:name" v="Lincoln"/>
    <tag k="addr:street:prefix" v="North"/>
    <tag k="addr:street:type" v="Avenue"/>
    <tag k="amenity" v="pharmacy"/>
    
      should be turned into:
    
    {...
    "address": {
        "housenumber": 5158,
        "street": "North Lincoln Avenue"
    }
    "amenity": "pharmacy",
    ...
    }
    
    - for "way" specifically:
    
      <nd ref="305896090"/>
      <nd ref="1719825889"/>
    
    should be turned into
    "node_refs": ["305896090", "1719825889"]
    """
    
    
    lower = re.compile(r'^([a-z]|_)*$')
    lower_colon = re.compile(r'^([a-z]|_)*:([a-z]|_)*$')
    problemchars = re.compile(r'[=\+/&<>;\'"\?%#$@\,\. \t\r\n]')
    addresschars = re.compile(r'addr:(\w+)')
    CREATED = [ "version", "changeset", "timestamp", "user", "uid"]
    OSM_FILE = 'jakarta_audit.osm'
    
    def shape_element(element):
        #node = defaultdict(set)
        node = {}
        if element.tag == "node" or element.tag == "way" :
            #create the dictionary based on exaclty the value in element attribute.
            node = {'created':{}, 'type':element.tag}
            for k in element.attrib:
                try:
                    v = element.attrib[k]
                except KeyError:
                    continue
                if k == 'lat' or k == 'lon':
                    continue
                if k in CREATED:
                    node['created'][k] = v
                else:
                    node[k] = v
            try:
                node['pos']=[float(element.attrib['lat']),float(element.attrib['lon'])]
            except KeyError:
                pass
            
            if 'address' not in node.keys():
                node['address'] = {}
            #Iterate the content of the tag
            for stag in element.iter('tag'):
                #Init the dictionry
    
                k = stag.attrib['k']
                v = stag.attrib['v']
                #Checking if indeed prefix with 'addr' and no ':' afterwards
                if k.startswith('addr:'):
                    if len(k.split(':')) == 2:
                        content = addresschars.search(k)
                        if content:
                            node['address'][content.group(1)] = v
                else:
                    node[k]=v
            if not node['address']:
                node.pop('address',None)
            #Special case when the tag == way,  scrap all the nd key
            if element.tag == "way":
                node['node_refs'] = []
                for nd in element.iter('nd'):
                    node['node_refs'].append(nd.attrib['ref'])
    #         if  'address' in node.keys():
    #             pprint.pprint(node['address'])
            return node
        else:
            return None
    
    
    def process_map(file_in, pretty = False):
        """
        Process the osm file to json file to be prepared for input file to monggo
        """
        file_out = "{0}.json".format(file_in)
        data = []
        with codecs.open(file_out, "w") as fo:
            for _, element in ET.iterparse(file_in):
                el = shape_element(element)
                if el:
                    data.append(el)
                    if pretty:
                        fo.write(json.dumps(el, indent=2)+"\n")
                    else:
                        fo.write(json.dumps(el) + "\n")
        return data
    
    def test():
    
        data = process_map(OSM_FILE)
        pprint.pprint(data[500])
    
    
    if __name__ == "__main__":
        test()