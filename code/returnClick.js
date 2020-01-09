$(document).keyup(function(event) {
    if ($("#gene_symbol").is(":focus") && (event.keyCode == 13)) {
        $("#go").click();
    }
});