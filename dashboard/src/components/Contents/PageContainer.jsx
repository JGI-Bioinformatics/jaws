import React from "react";

const ContentContainer = ({title, content}) => (
  <div className="MainPanel">
    <h1> {title} </h1>
    {content}
  </div>
)

export default ContentContainer;