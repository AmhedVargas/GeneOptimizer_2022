########Gene Optimizer 2022####
###Amhed Missael Vargas Velazquez
###avargas0lcg@gmail.com

###Short description:###
##By itself, the previous versions of the transgene builder could not allow multiple users due to the way that shiny is constructed, i.e. a single R session per app
##In theory, I could implement the package in a docker and mitigate the issues by creating a new session everytime a user connects it to. However, deploying a server only for that is a little cumbersome and not really required
##Similarly, I could use shiny proxy, but that also requires some configuration on the server's end and might obstruct the use we've implemented already.
##For now, the apparent most appropiate option is to use the package promises and start somewhat from scratch regarding the code. All right

#Load libraries
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(shinyjs)
library(DT)

####ExternalFunction
busyIndicator <- function(text = "Loading...", wait=2000) {
  shiny::tagList(
    shiny::div(class="loadbanner",id="loadmessage",text,img(src="elegans3.gif"))
    ,shiny::tags$script(sprintf(
      " setInterval(function(){
         if ($('html').hasClass('shiny-busy')) {
          setTimeout(function() {
            if ($('html').hasClass('shiny-busy')) {
              $('div.loadbanner').show()
            }
          }, %d)          
        } else {
          $('div.loadbanner').hide()
        }
      },100)
      ",wait)
    )
  ) 
}

# Define User interface
shinyUI(
    fluidPage(
      tags$head(
        tags$link(rel="stylesheet",type = "text/css", href="bootstrap.min.css")
      ),
      ###Loading message
      tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 60px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #D3D3D3;
               z-index: 105;
             }
          ")),
    ##Custom extra styles: single sliders background and title of navbar  
    tags$style(type = 'text/css', 
               ".js-irs-none .irs-single, .js-irs-none .irs-bar-edge, .js-irs-none .irs-bar {
                          background: transparent;
                          border-top-color: transparent;
                          border-bottom-color: transparent;
                          border-left-color: transparent;
                          border-right-color: transparent}
               .navbar-default .navbar-brand:hover {color: #ffffff;}
               
    .selectize-input {
        width: 600px;
        padding-top: 5px;
      }
               "),
    tags$head(tags$script(src ="sequence-viewer.bundle.js")),
    ###Extra scripts communicating java with shiny
    tags$script("
                var activitytimer;
                
                function checkActivity() {
  Shiny.setInputValue('clockactivity', Math.random());
  console.log(\" :stat\" + Math.random());
  activitytimer = setTimeout(function(){ checkActivity() }, 1000);
};

Shiny.addCustomMessageHandler('submitted-job', function(activity) {
console.log(activity + \" :stat1\");
        if(activity){
          checkActivity();
console.log(activity + \" :stat2\");
        }else{
console.log(activity + \" :statFInal\");
          clearTimeout(activitytimer);
        }

      });

Shiny.addCustomMessageHandler('JobStatusMessage', function(JobMessage) {
  document.getElementById(\"jobstat\").innerText = JobMessage;
      });
    "),
#Main tab pages
busyIndicator(),
    navbarPage(
      title=actionLink("link_to_tabpanel_sequenceadaptation", HTML("<b>Wormbuilder</b>")),
      windowTitle="WormBuilder transgenic tools",
        id = "panels",
        tabPanel("Sequence Adaptation",
                 mainPanel(
                     uiOutput("DynamicUserInterface")
                 )),
      ###About
      tabPanel("About",
               mainPanel(
                 h3("The app"),
                 HTML("<p align=\"justify\">Originally, the transgene builder app formed part of the <a href=\"https://wormbuilder.org/patc/\">PATC app</a> whose intention was to help
                      in the elaboration of transgenes less prone to germline silencing.
                      <br>
                      With this tool, we hope to leverage the time required to implement novel molecular tools yet to be seen in <i>C. elegans</i>.
                      </p>")
               )
      )
    ),
    HTML("<a href=\"https://syngenbio.kaust.edu.sa\">Syntetic genome biology laboratory @KAUST</a><br>"),
    HTML("<a href=\"http://www.wormbuilder.org/\">Wormbuilder</a><br>"),
    HTML("<a href=\"mailto:amhed.velazquez@kaust.edu.sa\">Contact us!</a>")
    
)
)

