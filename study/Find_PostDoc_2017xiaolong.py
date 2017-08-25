# -*- coding: utf-8 -*-
"""
Created on Fri May  5 15:15:08 2017

@author: k
"""
def getlist20170507():
    folder = 'C:\\Users\\k\\OneDrive\\mine\\future\\2017Find_PostDoc\\'
    cell = folder + 'Cell.xml'
    nature = folder + 'Nature.xml'
    science = folder + 'Science.xml'
    from lxml import etree
    
    
    CNS = list(etree.parse(cell).getroot()) + list(etree.parse(nature).getroot()) + list(etree.parse(science).getroot())
    
    #keep only those with abstract
    CNS_1 = [e for e in CNS if e.find('MedlineCitation/Article/Abstract') is not None]
    
    #keep only those author number is greater than min_authors
    min_authors = 4
    CNS_2 = [e for e in CNS_1 if len(e.findall('MedlineCitation/Article/AuthorList/Author'))>= min_authors]
    
    #the last author is considered as corresponding author, keep those with address in USA.
    country = 'USA'
    def authorfilter(paper, country = 'USA'):
        '''
        given an element of lxml.etree._Element of PubmedArticle
        return True if the corresponding author is in USA
        sometimes the author is not a human name, so there is no AffiliationInfo tag. skip that.
        if none of authors have AffiliationInfo, return false
        '''
        if len(paper.findall('MedlineCitation/Article/AuthorList/Author/AffiliationInfo/Affiliation')) == 0:
            return False
        authors = paper.findall('MedlineCitation/Article/AuthorList/Author')
        authors = [author for author in authors if len(author.findall('AffiliationInfo/Affiliation'))>0]
        correspond = authors[-1]
        Affiliation = correspond.findall('AffiliationInfo/Affiliation')
        if country in Affiliation[-1].text:
            return True
        return False
    CNS_3 = [e for e in CNS_2 if authorfilter(e,country)]
    
    #save PMID, PubDate, JournalName, title, Abstract, corresponding author name, info
    def getInfoFromPubmedArticle(paper):
        '''
        given an element of lxml.etree._Element of PubmedArticle
        return a tuple of PMID, PubDate, JournalName, title, Abstract, corresponding author name, email, address, author_number
        if country name is in the address of correspoinding author, total author numbers is greater or equal to min_authors, and there is an abstract, and at least one author have email address .
        return None if the standards are not satisfied.
        '''
        PMID = paper.find('MedlineCitation/PMID').text
        JournalName = paper.find('MedlineCitation/Article/Journal/Title').text
        PubDate = paper.find('MedlineCitation/Article/Journal/JournalIssue/PubDate')
        PubDate = '/'.join(e.text for e in PubDate)
        title = paper.find('MedlineCitation/Article/ArticleTitle').text
        Abstract = paper.find('MedlineCitation/Article/Abstract/AbstractText').text
        authors = paper.findall('MedlineCitation/Article/AuthorList/Author')
        authors = [author for author in authors if len(author.findall('AffiliationInfo/Affiliation'))>0]
        correspond = authors[-1]
        authorname = correspond.find('ForeName').text + ' ' + correspond.find('LastName').text
        address = correspond.findall('AffiliationInfo/Affiliation')[-1].text
        return PMID, JournalName, PubDate, title, Abstract, authorname, address
    
    paperSelected = [getInfoFromPubmedArticle(e) for e in CNS_3]
    #save paperSelected
    fout = open(folder + 'CNS_selected.txt','w',encoding='utf-8')
    for e in paperSelected:
        fout.write('\t'.join(e) +'\n')
    fout.close()
    
    