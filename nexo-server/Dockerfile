FROM node:7

RUN mkdir /nexo
WORKDIR /nexo

ADD . /nexo
RUN npm install

EXPOSE 3000

CMD ["node", "app.js"]
