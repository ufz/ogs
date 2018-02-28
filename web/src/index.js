import styles from './styles.css'

// os-selector
$("#btn-win").click(function(){
  $(".win").show();
  $(".linux").hide();
  $(".mac").hide();
  $("#btn-win").addClass("active");
  $("#btn-linux").removeClass("active");
  $("#btn-mac").removeClass("active");
});

$("#btn-linux").click(function(){
  $(".win").hide();
  $(".linux").show();
  $(".mac").hide();
  $("#btn-win").removeClass("active");
  $("#btn-linux").addClass("active");
  $("#btn-mac").removeClass("active");
});

$("#btn-mac").click(function(){
  $(".win").hide();
  $(".linux").hide();
  $(".mac").show();
  $("#btn-win").removeClass("active");
  $("#btn-linux").removeClass("active");
  $("#btn-mac").addClass("active");
});

$( document ).ready(function() {
  if ($(".win").length == 0) {
    $("#os-selector").hide();
  }
  else {
    $("#btn-win").click();
  }
});

// Site search
import algoliasearch from 'algoliasearch/lite';
var autocomplete = require('autocomplete.js/dist/autocomplete.jquery.js')

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
  // console.log(suggestion, dataset);
  window.location = "/" + suggestion.uri;
});
