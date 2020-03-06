$(document).keyup(function(event) {
    if ($("#gene_or_pathway_text").is(":focus") && (event.keyCode == 13)) {
        $("#gene_or_pathway_search").click();
    }
});