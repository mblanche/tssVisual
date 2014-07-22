setInterval(function(){
  if ($('html').attr('class')=='shiny-busy') {
    setTimeout(function() {
      if ($('html').attr('class')=='shiny-busy') {
        $('div.busy').show()
      }
    }, 1000)
  } else {
    $('div.busy').hide()
  }
}, 100)

// The UIs is not render until the data are loaded
// Use that to display a spinning disk
setInterval(function(){
    if (!$("#analysisType").is(":visible") ){
	setTimeout(function() {
	    if (!$("#analysisType").is(":visible") ){
		$('#dataLoader').show()
	    }
	}, 1000)
    } else {
	$('#dataLoader').hide()
    }
}, 100)
