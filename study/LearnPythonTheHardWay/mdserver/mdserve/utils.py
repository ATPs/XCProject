import re
import os
import markdown

def convert_markdown(doc):
    contents = open("docs/"+doc+".md").read()
    # convert it to markdown
    html = markdown.markdown(contents)
    # return it to the browser
    return html