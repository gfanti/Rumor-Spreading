# export graph

def export_gexf(filename,adjacency,source,infection_pattern,underlying_adjacency):
    # exports the graph described by adjacency as a GEXF file
    # Inputs:
    #   filename:   name of the file (no extension) to be saved
    #   adjacency:  the adjacency list for the graph
    #
    # Outputs:
    #   none, but creates a file named <filename>.gexf
    
    f = open('Gephi/'+filename+'.gexf','w')
    
    # write the header
    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write('<gexf xmlns="http://www.gexf.net/1.2draft"\n')
    f.write('    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n')
    f.write('    xsi:schemaLocation="http://www.gexf.net/1.2draft\n')
    f.write('          http://www.gexf.net/1.2draft/gexf.xsd"\n')
    f.write('    version="1.2">\n')
    
    # write the metadata
    f.write('    <meta lastmodifieddate="2009-03-20">\n')
    f.write('      <creator>Giulia</creator>\n')
    f.write('      <description>Adjacency graph for infection</description>\n')
    f.write('    </meta>\n')
    
    # build the graph
    f.write('    <graph defaultedgetype="undirected">\n')
    f.write('      <attributes class="node">\n')
    f.write('        <attribute id="0" title="Infection Status" type="string">\n')
    f.write('          <default>"Infected"</default>\n')
    f.write('        </attribute>\n')
    f.write('      </attributes>\n')
    
    
    # write the graph nodes
    f.write('      <nodes>\n')
    for i in range(len(adjacency)):
        if adjacency[i]:
            if i == source:
                f.write('        <node id="'+str(i)+'" label="'+str(i)+'">\n')
                f.write('          <attvalues>\n')
                f.write('            <attvalue for="0" value="Source"/> \n')
                f.write('          </attvalues> \n')
                f.write('        </node> \n')
            else:
                f.write('        <node id="'+str(i)+'" label="'+str(i)+'" />\n')
            if len(adjacency[i]) == 1:
                for j in underlying_adjacency[i]:
                    if not infection_pattern[j]:
                        # add an uninfected node
                        f.write('        <node id="'+str(j)+'" label="'+str(j)+'">\n')
                        f.write('          <attvalues>\n')
                        f.write('            <attvalue for="0" value="Uninfected"/> \n')
                        f.write('          </attvalues> \n')
                        f.write('        </node> \n')
                        
                        
    f.write('      </nodes>\n')
    
    # write the graph edges
    f.write('      <edges>\n')
    count = 0
    for i in range(len(adjacency)):
        if adjacency[i]:
            for target in adjacency[i]:
                if target > i:
                    f.write('        <edge id="'+str(count)+'" source="'+str(i)+'" target="'+str(target)+'" />\n')
                    count += 1
            if len(adjacency[i]) == 1:
                for target in underlying_adjacency[i]:
                    if not infection_pattern[target]:
                        # add an edge to an uninfected node
                        f.write('        <edge id="'+str(count)+'" source="'+str(i)+'" target="'+str(target)+'" />\n')
                        count += 1
    f.write('      </edges>\n')
    
    # end the graph
    f.write('    </graph>\n')
    
    # end the gexf file
    f.write('</gexf>\n')
    
    # close the file
    f.close()