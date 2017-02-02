import requests
import imp
import json
import os

# Load variables
if os.path.exists('./secret.py'):
    with open('./secret.py') as secretfile:
        vars = imp.load_source('vars', '', secretfile)
else:
    with open('./secret-from-env.py') as secretfile:
        vars = imp.load_source('vars', '', secretfile)


def send_request():
    # All news entries
    baseUrl = 'https://cdn.contentful.com/spaces'

    try:
        response = requests.get(
            url= baseUrl + '/' + vars.ogsSpace + '/entries',
            params={
                "content_type": "newsPost",
            },
            headers={
                "Authorization": "Bearer " + vars.accessToken,
            },
        )
        # print('Response HTTP Status Code: {status_code}'.format(
            # status_code=response.status_code))
        # print('Response HTTP Response Body: {content}'.format(
            # content=response.content))
        return response
    except requests.exceptions.RequestException:
        print('HTTP Request failed')


jsonData = send_request().json()

items = jsonData['items']

slugs = []
for item in items:
    id = item['sys']['id']
    slug = item['fields']['slug']
    title = item['fields']['title']
    slugs.append(slug)
    print id + " -> " + slug + ": " + title

with open('./../data/news.json', 'w') as outfile:
    json.dump(jsonData, outfile)

print slugs
with open('./../content/internal/news.md', 'w') as outfile:
    outfile.write("+++\nnews = [")
    for slug in slugs:
        outfile.write("\"%s\", " % slug)
    outfile.write("]\n+++\n")
