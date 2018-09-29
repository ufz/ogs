// Roughly based on https://github.com/algolia/autocomplete.js#jquery
// https://github.com/algolia/autocomplete.js/blob/master/LICENSE

var client = algoliasearch('4AHEU3QJQG', 'cda2754fe35733ffa31994a177725640')
var index = client.initIndex('docs.opengeosys.org');
$('#search-input').autocomplete({
    autoselect: true,
    hint: false,
    keyboardShortcuts: ['s']
}, [
  {
    source: $.fn.autocomplete.sources.hits(index, { hitsPerPage: 5 }),
    displayKey: 'title',
    templates: {
      suggestion: function(suggestion) {
        return suggestion._highlightResult.title.value;
      }
    }
  }
]).on('autocomplete:selected', function(event, suggestion, dataset) {
  window.location = "/" + suggestion.uri.toLowerCase();
});
